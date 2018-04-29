import csv
import logging
import os
from Bio import SeqIO
from ott_objects_defs.ploidbCommon import getPropertyFromFastaSeqHeader

__author__ = 'moshe'


"""
Handling multiple accessions.
Input:
	a fasta file of a single cluster
	BLAST result file
Output:
	The fasta file without multiple accessions

When a single species (taxon) has more than a single seq for this cluster, the script selects a single representative seq.
This is done by looking at the BLAST results and choosing the seq that its average distance from the rest in minimal
"""


logger = logging.getLogger('ploiDB-main')


# Given blast results the function will remove all rows where the alignment doesn't cover at least 50% of the sequences
def split_seqs_by_species(in_fasta_filename,work_dir):
	seq_data = SeqData()
	taxon_to_seqs_dict = dict()
	seqs = list(SeqIO.parse(in_fasta_filename, "fasta"))

	# Sorting the sequences according to the taxonid
	logger.debug("Splitting seqs by taxon - sorting %d taxons in %s into taxon files" % (len(seqs),in_fasta_filename))
	for seq_record in seqs:
		taxonid = getPropertyFromFastaSeqHeader(seq_record.description,"taxonid")
		gi = getPropertyFromFastaSeqHeader(seq_record.description,"gi")
		#seqid = getPropertyFromFastaSeqHeader(seq_record.description,"seqid") #Caused a bug in run Narcisuss when ITS cluster seqs headers were too long and so the seqid field got cut off.
		# Since it is the same now as the gi because it is actually the same value so:
		seqid = getPropertyFromFastaSeqHeader(seq_record.description,"gi")
		seq_data.blast_key_to_gi[seqid] = gi
		seq_data.blast_key_to_taxon[seqid] = taxonid
		if taxonid not in taxon_to_seqs_dict:
			taxon_to_seqs_dict[taxonid] = list()
		if taxonid not in seq_data.taxon_to_gis_dict:
			seq_data.taxon_to_gis_dict[taxonid] = list()
		taxon_to_seqs_dict[taxonid].append(seq_record)
		seq_data.taxon_to_gis_dict[taxonid].append(gi)

	logger.info("Writing each taxon to its own file. %d taxons were found" % (len(taxon_to_seqs_dict)))
	for taxonid in taxon_to_seqs_dict:
		taxon_filename = os.path.join(work_dir,taxonid + ".fasta")
		with open(taxon_filename,"w") as handle:
			SeqIO.write(taxon_to_seqs_dict[taxonid], handle, "fasta")
			seq_data.taxon_to_seqno[taxonid] = len(taxon_to_seqs_dict[taxonid])
			seq_data.taxon_to_filenames[taxonid] = taxon_filename

	logger.info("Done splitting by taxon")
	return seq_data


# Given blast results the function will remove all rows where the alignment doesn't cover at least 50% of the sequences
def split_blast_by_taxon_ids(in_blast_filename,work_dir,taxon_ids,seq_data):
	logger.info("Splitting BLAST results into separate files, per taxon")
	rows_dict = dict()

	for taxon_id in taxon_ids:
		rows_dict[taxon_id] = list()

	with open(in_blast_filename) as f:
		dr = csv.DictReader(f, delimiter='\t',fieldnames=['query_id','subject_id','pct_identity','align_len','mismatches',
														  'gap_openings','q_start','q_end','s_start','s_end','eval','bitscore'])

		logger.info("Filtering the blast results keeping only relevant results ofr multiple accessions")
		for next_row in dr:
			qid = next_row['query_id']
			sid = next_row['subject_id']
			taxon_qid = seq_data.blast_key_to_taxon[qid]
			gi_qid = seq_data.blast_key_to_gi[qid]
			taxon_sid = seq_data.blast_key_to_taxon[sid]
			gi_sid = seq_data.blast_key_to_gi[sid]

			if taxon_qid == taxon_sid and taxon_qid in taxon_ids:
				gis = seq_data.taxon_to_gis_dict[taxon_qid]
				if qid != sid and (gi_qid in gis or gi_sid in gis):
					rows_dict[taxon_qid].append(next_row)

		taxon_to_filenames = dict()

		logger.info("Writing the blast results into seperate files")
		for taxon_id in rows_dict:
			taxon_filename = os.path.join(work_dir,taxon_id + "-blast.csv")
			taxon_to_filenames[taxon_id] = taxon_filename
			with open(taxon_filename,"w", encoding='utf8', newline='') as handle:
				csv_writer = csv.writer(handle,delimiter='\t')
				for row in rows_dict[taxon_id]:
					row_to_write = [row['query_id'],row['subject_id'],row['pct_identity'],row['align_len'],row['mismatches'],
														  row['gap_openings'],row['q_start'],row['q_end'],row['s_start'],row['s_end'],row['eval'],row['bitscore']]
					csv_writer.writerow(row_to_write)

		logger.info("DONE Splitting BLAST results into separate files, per taxon")
		return taxon_to_filenames


def get_taxon_gis_key(gi1,gi2):
	key = "$".join([gi1,gi2])
	return key


def get_row_key(row):
	return "$".join([row['query_id'],row['subject_id']])


def get_relevant_seqids(fasta_filename):
	seqids = set()
	seqs = list(SeqIO.parse(fasta_filename, "fasta"))
	seqid_to_seq = dict()
	for seq in seqs:
		#seqid = getPropertyFromFastaSeqHeader(seq.description,"seqid")
		seqid = getPropertyFromFastaSeqHeader(seq.description,"gi")
		seqids.add(seqid)
		seqid_to_seq[seqid] = seq
	logger.info("Looking for GIs in %s - here is the GI list %s" % (fasta_filename,",".join(seqids)))

	return seqid_to_seq


def get_seq_with_max_average_blast_score(taxon_fasta_filename,taxon_blast_filename):
	seqids_to_seq = get_relevant_seqids(taxon_fasta_filename)

	logger.debug("Generating dictionary of bitscores between seqs according to %s" % taxon_blast_filename)
	with open(taxon_blast_filename) as f:
		dr = csv.DictReader(f, delimiter='\t',fieldnames=['query_id','subject_id','pct_identity','align_len','mismatches',
														  'gap_openings','q_start','q_end','s_start','s_end','eval','bitscore'])
		max_score_dict = dict()
		for row in dr:
			row_key = get_row_key(row)
			#logger.debug("Adding the following key %s" % row_key)
			if row_key not in max_score_dict:
				max_score_dict[row_key] = -1.0
			score = float(row['bitscore'])
			if max_score_dict[row_key] < score:
				max_score_dict[row_key] = score

	seqid_to_average_score = dict()
	missing_keys = list()
	for seqid in seqids_to_seq:
		average_bit_score = 0.0
		for other_seqid in seqids_to_seq:
			if seqid != other_seqid:
				key = get_taxon_gis_key(seqid,other_seqid)
				if key not in max_score_dict:
					missing_keys.append(key)
				else:
					average_bit_score = average_bit_score + max_score_dict[key]
		average_bit_score /= len(seqids_to_seq) - 1
		seqid_to_average_score[seqid] = average_bit_score
	if len(missing_keys) > 0:
		logger.error("Didn't find the following keys in blast file %s: %s" % (taxon_blast_filename, ",".join(missing_keys)))

	max_seqid = None
	max_average_bitscore = -1
	for seqid in seqid_to_average_score:
		# second check is done in order to make sure this method will always return the same seqid in case there are several seqs with the same average_bitscore
		if (max_average_bitscore < seqid_to_average_score[seqid]) or (max_average_bitscore == seqid_to_average_score[seqid] and seqid > max_seqid) :
			max_average_bitscore = seqid_to_average_score[seqid]
			max_seqid = seqid

	logger.info("Max average bitscore is %f for %s .Found the following average bit scores per GI %s" % (max_average_bitscore,max_seqid,seqid_to_average_score))
	return seqids_to_seq[max_seqid]

class SeqData:

	def __init__(self):
		self.taxon_to_seqno = dict()
		self.taxon_to_filenames = dict()
		self.taxon_to_gis_dict = dict()
		self.blast_key_to_gi = dict()
		self.blast_key_to_taxon = dict()

	def __str__(self):
		return "SeqData len(taxon_to_seqno)=%d len(self.taxon_to_filenames) =%d len(self.taxon_to_gis_dict)=%d" % (len(self.taxon_to_seqno),len(self.taxon_to_filenames),len(self.taxon_to_gis_dict))


def handle_multiple_accessions(in_seqs_filename,in_blast_filename,work_dir,out_fasta_filename):
	seq_data = split_seqs_by_species(in_seqs_filename,work_dir)
	logger.debug("Returned seq_data is %s" % seq_data)

	taxon_ids_for_multiple_acc = list()
	for taxon_id in seq_data.taxon_to_seqno:
		if seq_data.taxon_to_seqno[taxon_id] > 1:
			taxon_ids_for_multiple_acc.append(taxon_id)

	taxon_to_blast_filenames = split_blast_by_taxon_ids(in_blast_filename,work_dir,taxon_ids_for_multiple_acc,seq_data)

	logger.debug("About to iterate %d taxons, handling mult accession for each if needed" % (len(seq_data.taxon_to_filenames)))
	seqs_after_mult_acc = list()
	for taxon_id in seq_data.taxon_to_filenames:
		taxon_fasta_filename = seq_data.taxon_to_filenames[taxon_id]
		if taxon_id in taxon_ids_for_multiple_acc:
			logger.info("taxonid = %s - handling multiple accessions" % taxon_id)
			taxon_blast_filename = taxon_to_blast_filenames[taxon_id]
			representative_seq = get_seq_with_max_average_blast_score(taxon_fasta_filename,taxon_blast_filename)
		else:
			#logger.debug("taxonid = %s - no need for handling multiple accessions since there is just one seq" % taxon_id)
			with open(taxon_fasta_filename) as handle:
				temp_seqs = list(SeqIO.parse(handle,"fasta"))
				if len(temp_seqs) != 1:
					logger.error("Expected single seq but found %d for %s" % (len(temp_seqs), taxon_fasta_filename))
				representative_seq = temp_seqs[0]
		seqs_after_mult_acc.append(representative_seq)

	logger.info("Wrting %d seqs after handling multiple accessions to %s" % (len(seqs_after_mult_acc),out_fasta_filename))
	with open(out_fasta_filename, "w") as handle:
		SeqIO.write(seqs_after_mult_acc, handle, "fasta")
