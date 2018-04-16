from __future__ import division
import re
import os
from helper_functions import *
from pipeline_defs import *
from TreeConstruction import DistanceCalculator
from Bio import AlignIO


__author__ = 'Shiran'


def get_chars_content(sequence, sel_chars, all_chars):
	"""returns the fraction of sel_chars within all_chars in sequence"""
	#chars = a string of characters concatenated within square parenthesis to look for in sequence
	count = len(re.findall(sel_chars, sequence, re.I))
	total_base_count = len(re.findall(all_chars, sequence, re.I))
	fraction = count/total_base_count
	return fraction


def get_chars_statistic_per_alignment(msa_seqs, sel_chars, all_chars):
	seqs_chars_content = []

	for rec in msa_seqs:
		seq = rec._seq._data
		chars_content = get_chars_content(seq, sel_chars, all_chars)
		seqs_chars_content.append(chars_content)

	avg = get_avg(seqs_chars_content)
	std = get_std(seqs_chars_content)

	return avg, std


def remove_gaps_from_sequence(seq):
	return re.sub("[^a^g^c^t^A^G^C^T]+", "", seq, re.I)


def is_column_invariant(col_string):
	col_gapless = remove_gaps_from_sequence(col_string)
	#if something that is not gaps but not AGCT appears in the column
	if len(col_gapless) == 0:
		return 0
	#return true if the whole string is composed of repetitions of the first letter regardless of gaps
	return col_gapless == col_gapless[0]*len(col_gapless)


def get_p_invariant_sites(msa):
	msa_length = msa.get_alignment_length()
	invariant_sites = 0
	for col_i in range(0, msa_length):
		invariant_sites += is_column_invariant(msa.get_column(col_i))

	return invariant_sites/float(msa_length)


def calculate_column_entropy(col_string):
	# column_entropy = - sum(for every nucleotide x) {count(x)*log2(Prob(nuc x in col i))}
	col_gapless = remove_gaps_from_sequence(col_string).upper()
	col_entropy = 0
	for x in ['A', 'G', 'C', 'T']:
		count_x = str.count(col_gapless, x)
		if count_x == 0:
			entropy_x = 0
		else:
			prob_x = count_x/len(col_gapless)
			entropy_x = count_x*math.log2(prob_x)
		col_entropy += entropy_x

	return -col_entropy


def get_msa_avg_entropy(msa):
	msa_length = msa.get_alignment_length()
	sum_entropy = 0
	for col_i in range(0, msa_length):
		sum_entropy += calculate_column_entropy(msa.get_column(col_i))

	return sum_entropy/msa_length


def construct_distance_matrix(seqs_msa, output_path):
	dist_calc = DistanceCalculator()
	distance_matrix = dist_calc.get_distance(seqs_msa)

	names = distance_matrix.names
	dm_filepath = SEP.join([output_path, DM_PHYLIP_FILENAME])
	#don't overrun existing files
	if os.path.exists(dm_filepath):
		return dm_filepath
	fout = open(dm_filepath, 'w')
	fout.write('%d\n' % len(names))

	for name1 in names:
		fout.write(name1)
		for name2 in names:
			fout.write('\t%s' % distance_matrix[name1, name2])
		fout.write('\n')
	fout.close()

	return dm_filepath


def calculate_multinomial(msa):
	msa_length = msa.get_alignment_length()
	counts = []
	for col_i in range(0, msa_length):
		col = msa.get_column(col_i)
		cntA = col.count("A")
		cntC = col.count("C")
		cntG = col.count("G")
		cntT = col.count("T")
		counts.append((cntA, cntC, cntG, cntT))
	counts_set = set(counts)
	multinomial = 0
	for cnt in counts_set:
		c = counts.count(cnt)
		multinomial += c*math.log(c)
	multinomial -= msa_length*math.log(msa_length)
	return multinomial


def infer_pairwise_substitution_matrix(seq1, seq2):
	substitution_count_dictionary = {"AC": 0, "AG": 0, "AT": 0, "CG": 0, "CT": 0, "GT": 0, "AA": 0, "GG": 0, "CC": 0,
	                                 "TT": 0, "1s": 0, "2s": 0} #1s for one space vs nucleotide, 2s for 2 spaces
	pa_length = 0
	for i in range(0, len(seq1)):
		ch1 = min(seq1[i].upper(), seq2[i].upper())
		ch2 = max(seq1[i].upper(), seq2[i].upper())
		if ch1 not in ["A", "G", "C", "T"] and ch2 not in ["A", "G", "C", "T"]:
			substitution_count_dictionary["2s"] +=1
		elif ch1 in ["A", "G", "C", "T"] and ch2 in ["A", "G", "C", "T"]: #both nucleotides
			substitution_count_dictionary[ch1+ch2] += 1
			pa_length +=1
		else: #1 space
			substitution_count_dictionary["1s"] += 1

	return substitution_count_dictionary, pa_length


def compute_pairwise_substitution_rates(seq1, seq2):
	MATCH_SCORE = 1
	MISMATCH_SCORE = -1
	GAP_SCORE = 0

	substitution_count_dictionary, pa_length = infer_pairwise_substitution_matrix(seq1, seq2)
	if pa_length == 0:
		return 0, 0, 0
	transition_rate = float(substitution_count_dictionary["AG"] + substitution_count_dictionary["CT"]) / pa_length
	transversion_rate = float(substitution_count_dictionary["AC"] + substitution_count_dictionary["AT"] +
							  substitution_count_dictionary["CG"] + substitution_count_dictionary["GT"]) / pa_length
	match = (substitution_count_dictionary["AA"]+substitution_count_dictionary["CC"]+substitution_count_dictionary["GG"]+substitution_count_dictionary["TT"])*MATCH_SCORE
	mismatch = (substitution_count_dictionary["AC"]+substitution_count_dictionary["AG"]+substitution_count_dictionary["AT"] +
	           substitution_count_dictionary["CG"]+substitution_count_dictionary["CT"]+substitution_count_dictionary["GT"]+substitution_count_dictionary["1s"])*MISMATCH_SCORE
	gap = substitution_count_dictionary["1s"]*GAP_SCORE

	return transition_rate, transversion_rate, match+mismatch+gap


def calculate_avg_transition_transversion_rate(msa_seqs):
	transition_rates = []
	transversion_rates = []
	sop_score = 0
	for i in range(0, len(msa_seqs)-1):
		for j in range(i+1,len(msa_seqs)):
			seq1 = msa_seqs[i]._seq._data
			seq2 = msa_seqs[j]._seq._data
			transition_rate, transversion_rate, pair_sop = compute_pairwise_substitution_rates(seq1, seq2)
			transition_rates.append(transition_rate)
			transversion_rates.append(transversion_rate)
			sop_score += pair_sop

	return get_avg(transition_rates), get_avg(transversion_rates), sop_score


def get_msa_from_file(msa_file_path):
	#open file if exists
	if not os.path.exists(msa_file_path):
		return None
	try:
		msa = AlignIO.read(msa_file_path, PHYLIP_FORMAT)
	except:
		return None
	return msa


def get_seqs_values(msa_file_path, alignment_type):
	msa = get_msa_from_file(msa_file_path)
	if msa is None:
		return None
	seqs_msa = list(msa)
	msa_features = SequencesFeatures(len(msa), msa.get_alignment_length(), alignment_type)
	msa_features.gapless_matrix_size = sum([len(remove_gaps_from_sequence(str(seq.seq))) for seq in msa._records])
	msa_features.set_p_invariant_sites(get_p_invariant_sites(msa))
	msa_features.set_entropy(get_msa_avg_entropy(msa))
	msa_features.set_gc_values(*get_chars_statistic_per_alignment(seqs_msa, "[GC]", "[GCAT]"))
	msa_features.set_gap_values(*get_chars_statistic_per_alignment(seqs_msa, "[-?]", "[-?GCAT]"))

	freqs = []
	for nuc in ["A", "C", "G", "T"]:
		freq_nuc = get_chars_statistic_per_alignment(seqs_msa, nuc, "[GCAT]")[0]
		freqs.append(freq_nuc)
		msa_features.set_nucleotide_frequency(nuc, freq_nuc)

	msa_features.Echanges_per_site = 1/float(1-sum_of_squares(freqs))
	msa_features.freq_std = get_std(freqs)
	transition_avg, transversion_avg, sop_score = calculate_avg_transition_transversion_rate(seqs_msa)
	msa_features.set_substitution_rates(transition_avg, transversion_avg)
	msa_features.sop_score = sop_score

	return msa_features


#calculate_multinomial(get_msa_from_file("D:\\My Documents\\Desktop\\BIC\\new .phy"))