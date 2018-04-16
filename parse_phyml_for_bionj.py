from __future__ import division
from pipeline_defs import *

__author__ = 'Shiran'


def get_frequency(nuc1, nuc2, gtr_file_reader):
	freq = re.search("(?<=  "+nuc1+" <-> "+nuc2+") {2,4}[0-9\.]*", gtr_file_reader).group(0)
	freq = re.sub(" ", "", freq)
	return float(freq)

def read_matrix(gtr_file_reader):
	mat = re.search("\. Instantaneous rate matrix :\s+\[A\-+C\-+G\-+T\-+\](\s+.*[0-9]\.[0-9]{5} *){4}$", gtr_file_reader, re.M).group(0)
	matrix_reader = re.findall("\-?[0-9]\.[0-9]{5}", mat, re.M)
	matrix = [[], [], [], []]
	for i in range(0,16):
		element = matrix_reader[i]
		matrix[i//4].append(float(element))
	return matrix


def search_for_value(lookup_string, file_reader):
	value = re.search("(?<= " + lookup_string + ")[0-9\.]+", file_reader)
	return value.group(0)


def get_base_tree_features(hky_phyml_file_filepath, gtr_phyml_file_filepath):
	#open file if exists
	hky_file = open(hky_phyml_file_filepath, 'r')
	hky_file_reader = hky_file.read()

	gtr_file = open(gtr_phyml_file_filepath, 'r')
	gtr_file_reader = gtr_file.read()

	tree_features = BaseTreeFeatures()

	try:
		tree_features.parsimony = search_for_value("Parsimony: \t{4}", hky_file_reader)
		tree_features.tree_size = search_for_value("Tree size: \t{4}", hky_file_reader)
		tree_features.gamma = search_for_value("Gamma shape parameter: \t{2}", hky_file_reader)
		tree_features.prop_inv = search_for_value("Proportion of invariant: \t{2}", hky_file_reader)
		tree_features.Ts_ts_ratio = search_for_value("Transition/transversion ratio: \t{1}", hky_file_reader)

		freq_A = float(re.search("(?<=f\(A\)\= )[\.0-9]+",gtr_file_reader).group(0))
		freq_C = float(re.search("(?<=f\(C\)\= )[\.0-9]+",gtr_file_reader).group(0))
		freq_G = float(re.search("(?<=f\(G\)\= )[\.0-9]+",gtr_file_reader).group(0))
		freq_T = float(re.search("(?<=f\(T\)\= )[\.0-9]+",gtr_file_reader).group(0))
		nuc_freq = [freq_A, freq_C, freq_G, freq_T]

		tree_features.freq_AC = get_frequency("A", "C", gtr_file_reader)
		tree_features.freq_AG = get_frequency("A", "G", gtr_file_reader)
		tree_features.freq_AT = get_frequency("A", "T", gtr_file_reader)
		tree_features.freq_CG = get_frequency("C", "G", gtr_file_reader)
		tree_features.freq_CT = get_frequency("C", "T", gtr_file_reader)
		tree_features.freq_GT = get_frequency("G", "T", gtr_file_reader)
		tree_features.freqs_std = get_std([tree_features.freq_AC, tree_features.freq_AG, tree_features.freq_AT, tree_features.freq_CG, tree_features.freq_CT, tree_features.freq_GT])
		subs = [[None, tree_features.freq_AC, tree_features.freq_AG, tree_features.freq_AT],
				[None, None, tree_features.freq_CG, tree_features.freq_CT],
				[None, None, None, tree_features.freq_GT]]
		matrix = read_matrix(gtr_file_reader)

		i = 0
		j = 1
		while (matrix[i][j] == 0) and (i < 4):
			if j == 3:
				i += 1
				j = i+1
			else:
				j += 1

		tree_features.rate = matrix[i][j]/float(nuc_freq[j]*subs[i][j])
	except:
		return None
	finally:
		hky_file.close()
		gtr_file.close()

	return tree_features
