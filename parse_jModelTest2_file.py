__author__ = 'Shiran'

import re
import glob
#from pipeline_defs import *
import parse_newick_tree



TABLE_HEADERS = "Model\s+-lnL\s+K\s+TEST\s+delta\s+weight\s+cumWeight(\s+uDelta)?"
TREE_LINE = "Tree for the best TEST model = "



def get_test_features(test, file_reader):
	"""
	returns selected model, its score according to the current test table
	"""

	#create TestResults object
	test_features = TestResults(JMODELTEST,test)

	#extract tree features
	test_tree_line = re.sub("TEST", test, TREE_LINE)
	selected_model_tree = re.search("(?="+test_tree_line+").+", file_reader).group(0)
	test_tree_features = parse_newick_tree.parse_tree_string(selected_model_tree)

	#extract features from table
	test_table_headers = re.sub("TEST", test, TABLE_HEADERS)
	test_headers_match = re.search(test_table_headers + "\r?\n\-+\s\r?\n", file_reader, re.M)
	if test_headers_match is None:
		return None
	test_table_1st_line = re.search(".*", file_reader[test_headers_match.end():]).group(0)

	selected_model_line_features = re.split("\s+", test_table_1st_line)
	selected_model = selected_model_line_features[0]
	test_score = selected_model_line_features[3]

	selected_model = re.sub("\+", "", selected_model) #remove + from model name (.+I+G ==> .IG)

	test_features.set_values(selected_model, test_score, test_tree_features)
	return test_features


def get_jmt_values(genus_path, cluster_ID, tree_type_filename):
	#returns [AIC score, AIC model, AICc score, AICc model, BIC score, BIC model, taxa (count), alignment length]
	# from mrAIC output txt file

	#open file if exists
	jmt_file_path = glob.glob(SEP.join([genus_path, cluster_ID, "*" + tree_type_filename]))
	if len(jmt_file_path) == 0:
		return None
	jmt_file_path = jmt_file_path[0]
	jmt_file = open(jmt_file_path, 'r')
	file_reader = jmt_file.read()
	jmt_file.close()

	#extract ntaxa and length
	try:
		ntaxa = re.search("(?<=number of sequences: )[0-9]+", file_reader).group(0)
		nchars = re.search("(?<=number of sites: )[0-9]+", file_reader).group(0)
	except:
		return None
	try:
		jmt_aic_results = get_test_features(AIC, file_reader)
		jmt_aicc_results = get_test_features(AICc, file_reader)
		jmt_bic_results = get_test_features(BIC, file_reader)
	except:
		if re.search("Computation of likelihood scores discontinued \.{3}", file_reader):
			return None, None, None, ntaxa, nchars
		else:
			return -1

	return jmt_aic_results, jmt_aicc_results, jmt_bic_results, ntaxa, nchars