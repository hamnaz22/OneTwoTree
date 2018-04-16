import csv
import itertools
import logging
from Bio import SeqIO
from ploidbCommon import getPropertyFromFastaSeqHeader, write_standarize_fasta, exec_blast_all_vs_all

__author__ = 'moshe'

logger = logging.getLogger('ploiDB-main')


def blast_and_rev_fasta_seqs(original_fasta,blast_results_filename, reverse_seqs_fasta, short_desc_fasta,genus_name = None):
	exec_blast_all_vs_all(original_fasta,short_desc_fasta, blast_results_filename,genus_name)
	rev_fasta_seqs(original_fasta,blast_results_filename, reverse_seqs_fasta)


def rev_fasta_seqs(original_fasta,blast_results_filename, reverse_seqs_fasta):
	logger.info("About to reverse the fasta %s to %s using blast results %s" % (original_fasta, reverse_seqs_fasta,blast_results_filename))

	seqs = list(SeqIO.parse(original_fasta, "fasta"))
	rev_seqs = reverse_seqs(seqs, blast_results_filename)
	write_standarize_fasta(reverse_seqs_fasta,rev_seqs)


def reverse_seqs_in_fasta(cluster_context,genus_name, with_outgroup = True):
	if with_outgroup:
		original_fasta = cluster_context.all_seqs_with_outgroup_filename
	else:
		original_fasta = cluster_context.all_seqs_fasta_filename_no_multiple_accessions

	short_desc_fasta = cluster_context.all_seqs_short_desc
	blast_results_filename = cluster_context.blast_results_filename
	reverse_seqs_fasta = cluster_context.reverse_seqs_in_fasta

	blast_and_rev_fasta_seqs(original_fasta,blast_results_filename, reverse_seqs_fasta, short_desc_fasta, genus_name)


# Given a list of seqs, reverse all seqs if the direction is not the same for all of them.
# Direction of seqs is determined by BLAST all-vs-all results
def reverse_seqs(seqs,blast_results_filename):
	with open(blast_results_filename) as f:
		dr = csv.DictReader(f, delimiter='\t',fieldnames=['query_id','subject_id','pct_identity','align_len','mismatches',
														  'gap_openings','q_start','q_end','s_start','s_end','eval','bitscore'])
		groups = itertools.groupby(dr, lambda d: d['query_id'])
		max_group = None
		max_key = None
		for k,g in groups:
			g_as_list = list(g)
			if max_group is None:
				max_group = g_as_list
				max_key = k
			if len(max_group) < len(g_as_list):
				logger.debug("%d" % (len(max_group)))
				max_group = g_as_list
				max_key = k

		logger.debug("largest group is %s size in blast results is %d" % (max_key,len(max_group)))

		# TODO: consider something like this:
		subjects = itertools.groupby(max_group, lambda d: d['subject_id'])
		seqs_to_reverse = set()
		for subject,subject_rows in subjects:
			logger.debug("checking subject %s" % subject)
			subject_rows_list = list(subject_rows)
			max_length_row = max(subject_rows_list, key=lambda d: int(d['align_len']))
			logger.debug("Max row out of %d has length of %s (query = %s)" % (len(subject_rows_list),max_length_row['align_len'],max_length_row['query_id']))
			start_end_diff = (float(max_length_row['s_end']) - float(max_length_row['s_start']))
			is_reversed = start_end_diff < 0
			if is_reversed:
				logger.debug("For %s %s the diff is %d and will be reversed" % (max_length_row['query_id'],max_length_row['subject_id'],start_end_diff))
				seqs_to_reverse.add(max_length_row['subject_id'])

	logger.info("The following seqs will be reversed %s" % ",".join(seqs_to_reverse))
	reversed_seqs = list()
	reversed_seqids = list()
	for seq in seqs:
		#logger.debug("type is %s" % type(seq))
		#seqid = getPropertyFromFastaSeqHeader(seq.description,"seqid")
		seqid = getPropertyFromFastaSeqHeader(seq.description,"gi")
		if seqid in seqs_to_reverse:
			rev_seq = seq.reverse_complement()
			rev_seq.description = seq.description
			reversed_seqs.append(rev_seq)
			reversed_seqids.append(seqid)
		else:
			reversed_seqs.append(seq)

	logger.info("The following seqs were actually reversed %s" % ",".join(reversed_seqids))

	return reversed_seqs
