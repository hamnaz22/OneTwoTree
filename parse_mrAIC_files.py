__author__ = 'Shiran'

import re
import glob
import parse_newick_tree
from pipeline_defs import *


def extract_model_scores(test_criteria, file_reader):
	values_cursor = re.search("Model\tdf\tlnL\t\t" + test_criteria + "\t", file_reader, re.MULTILINE)
	if values_cursor is None and test_criteria == "AICc":
		return None
	model_optimal_line = re.search("(?<=\n).+(?=\n)", file_reader[values_cursor.end(0):]).group(0)
	model_score = str.split(model_optimal_line, "\t")[3]

	return model_score


def extract_model_parameters(test_criteria, file_reader):
	model_description_block = re.search("\[Mrbayes block for the best " + test_criteria + " model.*?\nEND;",
										file_reader, re.S)
	if model_description_block is None and test_criteria == "AICc":
		return [None]*3
	model_description = model_description_block.group(0)
	nst_val = re.search("(?<=nst=)[126]",model_description, re.M).group(0)
	rates_val = re.search("(?<=rates=)[a-z]*(?=;\n)",model_description, re.M).group(0)

	#infer I and gamma from rates_val
	if rates_val == "equal":
		is_inv = False
		is_gamma = False
	elif rates_val == "propinv":
		is_inv = True
		is_gamma = False
	elif rates_val == "gamma":
		is_inv = False
		is_gamma = True
	elif rates_val == "invgamma":
		is_inv = True
		is_gamma = True

	return [is_inv, is_gamma]


def get_mraic_values(genus_path, cluster_ID):
	#returns [AIC score, AIC model, AICc score, AICc model, BIC score, BIC model, taxa (count), alignment length]
	# from mraic output txt file

	#open file if exists
	mraic_file_path = glob.glob(SEP.join([genus_path, cluster_ID, "*" + MRAIC_TXT_SUFFIX]))
	if len(mraic_file_path) == 0:
		return None
	mraic_file_path = mraic_file_path[0]
	mraic_file = open(mraic_file_path, 'r')
	file_reader = mraic_file.read()

	#extract taxa and length
	taxa_length_line = re.search("Input data from file.*\n", file_reader).group(0)
	ntaxa = re.search("[0-9]+(?= taxa)", taxa_length_line).group(0)
	nchars = re.search("[0-9]+(?= characters)", taxa_length_line).group(0)

	#extract selected models
	aic_model = re.search("(?<=Minimum AIC  model: ).+", file_reader).group(0)
	aicc_model = re.search("(?<=Minimum AICc model: ).+", file_reader).group(0)
	bic_model = re.search("(?<=Minimum BIC  model: ).+", file_reader).group(0)

	#extract model scores
	aic_score = extract_model_scores("AIC", file_reader)
	aicc_score = extract_model_scores("AICc", file_reader)
	if aicc_score is None:
		aicc_score = aic_score
	bic_score = extract_model_scores("BIC", file_reader)

	#extract model parameters
	#aic_parameters = extract_model_parameters("AIC", file_reader)
	#aicc_parameters = extract_model_parameters("AICc", file_reader)
	#bic_parameters = extract_model_parameters("BIC", file_reader)

	aic_tree_properties = parse_newick_tree.parse_mraic_tree_file(SEP.join([genus_path, cluster_ID]), "AIC")
	aicc_tree_properties = parse_newick_tree.parse_mraic_tree_file(SEP.join([genus_path, cluster_ID]), "AICc")
	bic_tree_properties = parse_newick_tree.parse_mraic_tree_file(SEP.join([genus_path, cluster_ID]), "BIC")

	mraic_file.close()

	#set in test_result struct
	mraic_aic_results = TestResults(MR_AIC, AIC)
	mraic_aic_results.set_values(aic_model, aic_score, aic_tree_properties)
	mraic_aicc_results = TestResults(MR_AIC, AICc)
	mraic_aicc_results.set_values(aicc_model, aicc_score, aicc_tree_properties)
	mraic_bic_results = TestResults(MR_AIC, BIC)
	mraic_bic_results.set_values(bic_model, bic_score, bic_tree_properties)

	return mraic_aic_results, mraic_aicc_results, mraic_bic_results, ntaxa, nchars
