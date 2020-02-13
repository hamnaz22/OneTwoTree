import argparse
import csv
import logging
import sys
from Bio import SeqIO
import itertools

################################################################################################################################################
# The purpose of this script is that given a BLAST all-vs-all file and the fasta file used for BLAST, it outputs a new BLAST resutls file
# after filtering short sequences.
# This is useful for cases where BLAST of 2 sequences gave a very good evalue/bitscore but actually one of the sequences is much shorted than
# the other (e.g. seq1 is contained in seq2). So even though they have good evalue we want to ignore them in our results since most likely
# they shouldn't be clustered together
################################################################################################################################################
__author__ = 'moshe'

ALIGN_LEN_THRESHOLD = 0.35
SQ_LEN_MIM_RATIO = 0.5
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


def noDuple_Align_length_calc(EndPoint_StartKey_dict):

	calc_length=0;indx=0
	#logger.debug("Calc real Align length according to dict:")
	#logger.debug(EndPoint_StartKey_dict)
	sorted_start_points = list(EndPoint_StartKey_dict.keys())
	#logger.debug(sorted_start_points)
	indx_max = len(sorted_start_points)
	if indx_max == 1:
		start_point = sorted_start_points[0]
		calc_length+=abs(EndPoint_StartKey_dict[start_point]-start_point)+1
		#logger.debug("return the follwoing align length:")
		#logger.debug(calc_length)
		return calc_length
	sorted_start_points = sorted(sorted_start_points)
	#logger.debug("sorted_start_points:")
	#logger.debug(sorted_start_points)
	#indx_max = len(EndPoint_StartKey_dict) #Number of pairs to check
	#calc first align length:
	start_point = sorted_start_points[0]
	calc_length+=abs(EndPoint_StartKey_dict[start_point]-start_point)+1
	indx+=1
	while indx < indx_max:
		cur_start_point = sorted_start_points[indx]
		prev_end_point = EndPoint_StartKey_dict[sorted_start_points[indx-1]]
		if prev_end_point >= cur_start_point:
			cur_contribution = EndPoint_StartKey_dict[sorted_start_points[indx]] - prev_end_point
		else:
			cur_contribution = EndPoint_StartKey_dict[sorted_start_points[indx]] - cur_start_point + 1
		#logger.debug("Adding contribution: %s" %str(cur_contribution))
		calc_length+=cur_contribution
		indx+=1
	#logger.debug("return the follwoing align length:")
	#logger.debug(calc_length)
	return calc_length


# Given blast results the function will remove all rows where the alignment doesn't cover at least 50% of the sequences
def filter_seqs_by_match_length(in_blast_filename,out_blast_filename,in_fasta_filename,filet_ratio_threshold):

	seqs_dict = SeqIO.index(in_fasta_filename, "fasta")

	with open(in_blast_filename) as f,open(out_blast_filename, mode="w", encoding='utf8', newline='') as out_handle:
		dr = csv.DictReader(f, delimiter='\t',fieldnames=['query_id','subject_id','pct_identity','align_len','mismatches',
														  'gap_openings','q_start','q_end','s_start','s_end','eval','bitscore'])

		# Grouping by queryid and subjectid
		groups = itertools.groupby(dr, lambda d: "$".join(sorted(([d['query_id'],d['subject_id']]))))

		csv_writer = csv.writer(out_handle, delimiter='\t')

		for key,group in groups:
			group_as_list = list(group)
			qkey = group_as_list[0]['query_id']
			skey = group_as_list[0]['subject_id']
			qlen = len(seqs_dict[qkey])
			slen = len(seqs_dict[skey])

			total_align_len_for_pair = 0
			query_ranges_dict ={}
			src_ranges_dict = {}
			for row in group_as_list:
				align_len = float(row['align_len'])
				#Collect all relevant data for alignment length:
				query_start = int(row['q_start'])
				query_end = int(row['q_end'])
				if query_start > query_end:
					query_start,query_end = query_end,query_start
					#temp_var = query_start
					#query_start = query_end
					#query_end = temp_var
				src_start = int(row['s_start'])
				src_end = int(row['s_end'])
				if src_start > src_end:
					src_start,src_end = src_end,src_start
				query_ranges_dict[query_start]=query_end
				src_ranges_dict[src_start]=src_end

				total_align_len_for_pair += align_len
				# Just for verifying that all the subject and query are the same in this group
				if qkey != row['query_id'] or skey != row['subject_id']:
					logger.critical("Mismatch between expected subject/query %s %s" % (qkey,skey))

			logger.debug("src and query alignment ranges dict: quer %s, src %s" %(qkey,skey))
			logger.debug(query_ranges_dict)
			logger.debug(src_ranges_dict)
			total_align_len_for_q_pair = noDuple_Align_length_calc(query_ranges_dict)
			total_align_len_for_s_pair = noDuple_Align_length_calc(src_ranges_dict)
			# Check all results for overlap and remove redundant:
			#x = range(1,10)
			#y = range(8,20)
			#xs = set(x)
			#xs.intersection(y)


			qlen_total_align = (total_align_len_for_q_pair/qlen)
			slen_total_align = (total_align_len_for_s_pair/slen)


			# Removing all rows of subject and query in which the alignment doesn't cover at least 50%
			#if (slen / qlen) < SQ_LEN_MIM_RATIO or (qlen / slen) < SQ_LEN_MIM_RATIO:
			#	logger.info("Skipping rows with qid = %s and sid = %s. lengths are qlen=%f slen=%f and ratios are =%f and %f" %
			#				 (qkey,skey,qlen,slen,slen / qlen,qlen / slen))
			if qlen_total_align < filet_ratio_threshold or slen_total_align < filet_ratio_threshold:
				logger.info("Skipping rows with qid = %s and sid = %s. lengths are qlen=%f slen=%f total_align_len_for_pair=%f qlen_total_align=%f and slen_total_align=%f (filter_ratio=%f)" %
							 (qkey,skey,qlen,slen,total_align_len_for_pair,qlen_total_align,slen_total_align,filet_ratio_threshold))
			else:
				for row in group_as_list:
					csv_writer.writerow([row['query_id'],row['subject_id'],row['pct_identity'],row['align_len'],row['mismatches'],
										row['gap_openings'],row['q_start'],row['q_end'],row['s_start'],row['s_end'],row['eval'],row['bitscore']])


#			if qlen_total_align < ALIGN_LEN_THRESHOLD or slen_total_align < ALIGN_LEN_THRESHOLD:
#				logger.info("Skipping rows with qid = %s and sid = %s. lengths are qlen=%f slen=%f total_align_len_for_pair=%f qlen_total_align=%f and slen_total_align=%f" %
#							 (qkey,skey,qlen,slen,total_align_len_for_pair,qlen_total_align,slen_total_align))
#			elif (slen / qlen) < SQ_LEN_MIM_RATIO or (qlen / slen) < SQ_LEN_MIM_RATIO:
#				logger.info("Skipping rows with qid = %s and sid = %s. lengths are qlen=%f slen=%f and ratios are =%f and %f" %
#							 (qkey,skey,qlen,slen,slen / qlen,qlen / slen))
#			else:
#				for row in group_as_list:
#					csv_writer.writerow([row['query_id'],row['subject_id'],row['pct_identity'],row['align_len'],row['mismatches'],
#										row['gap_openings'],row['q_start'],row['q_end'],row['s_start'],row['s_end'],row['eval'],row['bitscore']])



if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='The scripts takes a FASTA file with similar seqs (according to BLAST), and reverse sequence if needed')
	parser.add_argument('--in-filename','-i', help='Input blast results', required=True)
	parser.add_argument('--in-fasta-filename','-f', help='Input FASTA with seqs', required=True)
	parser.add_argument('--out-filename','-o', help='Output blast results', required=True)
	parser.add_argument('--filter-blast-ratio','-fbr', help='Filter ratio - seq vs total alignment', required=True)

	# -fbr $filterBlastRatio

	args = parser.parse_args()

	initCommandLineLogger()
	logger.info("About to remove short/long seqs from fasta %s to %s (filet ratio threshold = %s)" % (args.in_filename,args.out_filename,args.filter_blast_ratio))
	filter_seqs_by_match_length(args.in_filename,args.out_filename, args.in_fasta_filename, float(args.filter_blast_ratio))
