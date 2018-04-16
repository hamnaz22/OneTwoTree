import argparse
import csv
import logging
import sys
from Bio import SeqIO
import itertools
from collections import Counter
################################################################################################################################################
# The purpose of this script is that given a BLAST all-vs-all file and the fasta file used for BLAST, it outputs a new BLAST resutls file
# after filtering short sequences.
# This is useful for cases where BLAST of 2 sequences gave a very good evalue/bitscore but actually one of the sequences is much shorted than
# the other (e.g. seq1 is contained in seq2). So even though they have good evalue we want to ignore them in our results since most likely
# they shouldn't be clustered together
################################################################################################################################################
__author__ = 'moshe'

ALIGN_LEN_THRESHOLD = 0.5
SQ_LEN_MIM_RATIO = 0.5
logger = logging.getLogger('ploiDB-main')

class MultiIndexDict:
    def __init__(self, *indexes):
        self._indexes = indexes
    def __getitem__(self, key):
        for idx in self._indexes:
            try:
                return idx[key]
            except KeyError:
                pass
        raise KeyError("{0} not found".format(key))

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


# Given blast results the function will remove all rows where the alignment doesn't cover at least 50% of the sequences
def filter_seqs_by_match_length(in_blast_filename,out_blast_filename,in_fasta_filename,in_fasta_seq_filename):

	seqs_dict_s = SeqIO.index(in_fasta_filename, "fasta")
	seqs_dict_q = SeqIO.index(in_fasta_seq_filename, "fasta")
	for key in seqs_dict_s.keys():
		key_text = key
		break
	print(key_text)
	#add_txt = '|taxonid' + key_text.split('|taxonid')[1]
	#print(add_txt)

	with open(in_blast_filename) as f,open(out_blast_filename, mode="w", encoding='utf8', newline='') as out_handle:
		dr = csv.DictReader(f, delimiter='\t',fieldnames=['query_id','subject_id','pct_identity','align_len','mismatches',
														  'gap_openings','q_start','q_end','s_start','s_end','eval','bitscore'])

		# Grouping by queryid and subjectid
		#groups = itertools.groupby(dr, lambda d: "$".join(sorted(([d['query_id'],d['subject_id']]))))
		groups = itertools.groupby(dr, lambda d: "$".join(sorted((d['subject_id']))))

		csv_writer = csv.writer(out_handle, delimiter='\t')

		for key,group in groups:
			group_as_list = list(group)
			qkey = group_as_list[0]['query_id']
			skey = group_as_list[0]['subject_id'] #+ add_txt
			#qlen = len(seqs_dict[qkey])
			#slen = len(seqs_dict[skey])
			qlen = len(seqs_dict_q[qkey])
			slen = len(seqs_dict_s[skey])

			total_align_len_for_pair = 0
			for row in group_as_list:
				align_len = float(row['align_len'])
				total_align_len_for_pair += align_len
				# Just for verifying that all the subject and query are the same in this group
				#if qkey != row['query_id'] or skey.replace(add_txt,'') != row['subject_id']:
				if qkey != row['query_id'] or skey != row['subject_id']:
					logger.critical("Mismatch between expected subject/query %s %s" % (qkey,skey))

			qlen_total_align = (total_align_len_for_pair/qlen)
			slen_total_align = (total_align_len_for_pair/slen)


			# Removing all rows of subject and query in which the alignment doesn't cover at least 50%
			if qlen_total_align < ALIGN_LEN_THRESHOLD or slen_total_align < ALIGN_LEN_THRESHOLD:
				logger.info("Skipping rows with qid = %s and sid = %s. lengths are qlen=%f slen=%f total_align_len_for_pair=%f qlen_total_align=%f and slen_total_align=%f" %
							 (qkey,skey,qlen,slen,total_align_len_for_pair,qlen_total_align,slen_total_align))
			elif (slen / qlen) < SQ_LEN_MIM_RATIO or (qlen / slen) < SQ_LEN_MIM_RATIO:
				logger.info("Skipping rows with qid = %s and sid = %s. lengths are qlen=%f slen=%f and ratios are =%f and %f" %
							 (qkey,skey,qlen,slen,slen / qlen,qlen / slen))
			else:
				for row in group_as_list:
					csv_writer.writerow([row['query_id'],row['subject_id'],row['pct_identity'],row['align_len'],row['mismatches'],
										row['gap_openings'],row['q_start'],row['q_end'],row['s_start'],row['s_end'],row['eval'],row['bitscore']])



if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='The scripts takes a FASTA file with similar seqs (according to BLAST), and reverse sequence if needed')
	parser.add_argument('--in-filename','-i', help='Input blast results', required=True)
	parser.add_argument('--in-fasta-seq-filename','-i_seq', help='Input FASTA with seq', required=True)
	parser.add_argument('--in-fasta-db-filename','-f', help='Input FASTA with database seqs', required=True)
	parser.add_argument('--out-filename','-o', help='Output blast results', required=True)

	args = parser.parse_args()

	initCommandLineLogger()
	logger.info("About to remove short/long seqs from fasta %s to %s" % (args.in_filename,args.out_filename))
	filter_seqs_by_match_length(args.in_filename,args.out_filename, args.in_fasta_db_filename,args.in_fasta_seq_filename)


