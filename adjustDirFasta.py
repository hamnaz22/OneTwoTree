import argparse
import logging
import shutil
from Bio import SeqIO, Phylo
import os
import sys
from ploidbCommon import getPropertyFromFastaSeqHeader, exec_external_command_redirect_output, write_standarize_fasta

"""
 Given a list of fasta files, reverse all seqs in them if the direction is not the same for all seqs.
 This is done using mafft  --adjustdirection option.

Input: a list of fasta files (seperated by commas when calling from command line)
output: output is written to input files. The original files are saved with a ".bak" extenstion

 Given a fasta with multiple seqs, the script
 (1) executes mafft with --adjustdirection  option
 (2) examines mafft output in order to find seqs which their header was added a "_R" and if so, reverse them in the original fasta file
"""

__author__ = 'moshe'

logger = logging.getLogger('ploiDB-main')

def initCommandLineLogger():
	logger.setLevel(logging.DEBUG)
	# create console handler and set level to debug
	ch = logging.StreamHandler(sys.stdout)
	ch.setLevel(logging.DEBUG)
	# create formatter
	formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
	# add formatter to ch
	ch.setFormatter(formatter)
	# add ch to logger
	logger.addHandler(ch)


def adjust_direction_fasta_using_mafft(fasta_list, backup_ext = ".bak"):
	initial_fasta = fasta_list[0]
	work_dir, tail = os.path.split(initial_fasta)

	initial_fasta_bak = initial_fasta + backup_ext
	fasta_aligned_temp = os.path.join(work_dir,"adjustdir.tmp")
	fasta_aligned = os.path.join(work_dir,"adjustdir.aln")
	shutil.copy(initial_fasta,initial_fasta_bak)
	logger.debug("Backing %s to %s" % (initial_fasta,initial_fasta_bak))

	initial_fasta_aln_out = initial_fasta + ".aln.out"
	initial_fasta_aln_err = initial_fasta + ".aln.err"
	mafft_command = "mafft --adjustdirection  --maxiterate 500 --localpair %s > %s" % (initial_fasta, fasta_aligned_temp)
	if(exec_external_command_redirect_output(mafft_command,initial_fasta_aln_out,initial_fasta_aln_err)) != 0:
		logger.debug("ERROR: mafft_command Failed in adjust_direction_fasta_using_mafft")
		return
	#exec_external_command_redirect_output(mafft_command,initial_fasta_aln_out,initial_fasta_aln_err)
	logger.debug("Backing %s to %s" % (fasta_aligned_temp,fasta_aligned))
	shutil.copy(fasta_aligned_temp,fasta_aligned)
	for next_fasta in fasta_list[1:]:
		logger.debug("Processing %s" % next_fasta)
		next_fasta_bak = next_fasta + backup_ext
		next_fasta_aln_out = next_fasta + ".aln.out"
		next_fasta_aln_err = next_fasta + ".aln.err"
		logger.debug("Backing %s to %s" % (next_fasta,next_fasta_bak))
		shutil.copy(next_fasta,next_fasta_bak)
		mafft_command = "mafft --adjustdirection --addfragments %s --multipair %s > %s" % (next_fasta,fasta_aligned,fasta_aligned_temp)
		exec_external_command_redirect_output(mafft_command,next_fasta_aln_out,next_fasta_aln_err)
		shutil.copy(fasta_aligned_temp,fasta_aligned)

	PREFIX_FOR_REV_SEQ = "_R_gi"
	seqs = list(SeqIO.parse(fasta_aligned,"fasta"))
	gi_to_rev = list()
	for seq in seqs:
		if seq.description.startswith(PREFIX_FOR_REV_SEQ):
			new_desc = seq.description.replace(PREFIX_FOR_REV_SEQ, "gi")
			gi = getPropertyFromFastaSeqHeader(new_desc, "gi")
			gi_to_rev.append(gi)
			logger.debug("GI to rev %s" % gi)

	for next_fasta in fasta_list:
		logger.debug("Reversing seqs in %s" % next_fasta)
		seqs = list(SeqIO.parse(next_fasta, "fasta"))
		reversed_seqs = list()
		for seq in seqs:
			gi = getPropertyFromFastaSeqHeader(seq.description, "gi")
			if gi in gi_to_rev:
				rev_seq = seq.reverse_complement()
				rev_seq.description = seq.description
				reversed_seqs.append(rev_seq)
				logger.debug("Seq %s will be reversed" % rev_seq.description)
			else:
				reversed_seqs.append(seq)

		write_standarize_fasta(next_fasta,reversed_seqs)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='The scripts takes a FASTA file with similar seqs (according to BLAST), and reverse sequence if needed')
	parser.add_argument('--fasta-in-filenames','-i', help='Input FASTA with seqs', required=True)

	args = parser.parse_args()

	initCommandLineLogger()

	fastas_list = args.fasta_in_filenames.split(",")
	adjust_direction_fasta_using_mafft(fastas_list)


