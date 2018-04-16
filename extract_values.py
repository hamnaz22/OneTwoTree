import os
import csv
import parse_mrAIC_files
import parse_jModelTest2_file
import parse_MSA_file
import parse_newick_tree
from pipeline_defs import *
from helper_functions import *
import parse_phyml_for_bionj
import logging
import argparse

__author__ = 'Shiran'

"""
collects data per genus-cluster from mrAIC input and output files and fasta alignment files
"""


def tree_categories_with_source(source):
	tree_categories = ["total tree distance", "avg height", "height std", "max height" ]
	for i in range(len(tree_categories)):
		tree_categories[i] = " ".join([source, tree_categories[i]])
	return tree_categories

def categories_with_source(source):
	categories = ["model", "base-model", "score", "Inv", "Gamma", "F", "substitution count"]
	for i in range(len(categories)):
		categories[i] = " ".join([source, categories[i]])
	return categories + tree_categories_with_source(source)

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Runs mrAIC per genus')
	parser.add_argument('--db_type', '-db', help='Which DB- trpipe or AA_DB', required=True, default=TRPIPE_DB_NAME)
	parser.add_argument('--output_file', '-o', help='Output file name', required=True)
	args = parser.parse_args()

	if args.db_type == TRPIPE_DB_NAME:
		source_path = TRPIPE_OUTPUT_DIRECTORY
		alignment_type = "nucleotide"
		source_dict = import_csv_to_dict(TRPIPE_CLUSTER_INFO)
	elif args.db_type == AA_DB_NAME:
		source_path = AA_OUTPUT_DIRECTORY
		alignment_type = "proteinStruc"
	elif args.db_type == PANDIT_DB_NAME:
		source_path = PANDIT_OUTPUT_DIRECTORY
		alignment_type = "codon"
	elif args.db_type == SELETOME_DB_NAME:
		source_path = SELECTOME_OUTPUT_DIRECTORY
		alignment_type = "codon"
	elif args.db_type == RNA_DB_NAME:
		source_path = RNA_OUTPUT_DIRECTORY
		alignment_type = "secRNA"
	dest_path = SEP.join([FEATURES_SUMMARY_OUTPUT_PATH, args.output_file])

	dataFile = open(dest_path, 'w', newline='')
	csv_writer = csv.writer(dataFile, dialect='excel')

	##prepare an instance of a cluster features object for csv table headers
	#D = PhyloObject("instance","instance").convert_to_dictionary()
	#
	#dataFile = open(dest_path, 'w', newline='')
	##dataFile.write(u'\ufeff'.encode('utf8'))
	#csv_writer = csv.DictWriter(dataFile, D.keys())
	#csv_writer.writeheader()

	dAIc_categories = ["dAICc0", "dAICc1","dAICc2","dAICc3","dAICc4","dAICc5","dAICc6","dAICc7","dAICc8","dAICc9","dAICc10","dAICc11"]
	ddAIc_categories = ['p11p10', 'p9p4', 'p10p9', 'p9p5', 'p10p7', 'p10p6', 'p10p5']
	csv_writer.writerow(["ALN_ID"]
						+ categories_with_source("jMT_BioNJ-AIC") + categories_with_source("jMT_BioNJ-AICc") + categories_with_source("jMT_BioNJ-BIC")
						+ categories_with_source("jMT_ML-AIC") + categories_with_source("jMT_ML-AICc") + categories_with_source("jMT_ML-BIC") +
						 ["alignment_type", "taxa", "length", "klogn", "gapless_matrix_size", "GC-content_avg(per)", "GC-content_std(per)", "gap-content_avg(per)",
						 "gap-content_std(per)", "aln_transition_avg", "freq_A", "freq_C", "freq_G", "freq_T", "freq_std", "Exp_changes_per_site", "inv_sites", "aln_entropy", "sop_score"] +
						 dAIc_categories + ddAIc_categories + ["parsimony", "gamma", "p-Inv",
						 "Tstv", "sub_AC/GT", "sub_AG/GT", "sub_AT/GT", "sub_CG/GT", "sub_CT/GT", "substitution_std",
						 "mean_instantaneous_rate"] + tree_categories_with_source("PhyML_HKY") + ["source", "type"])

	logger = logging.getLogger('model-selection-features')
	init_commandline_logger(logger)
	directories = os.listdir(source_path)

	for genus in directories:
		if os.path.isfile(SEP.join([source_path, genus])):
			continue
		genus_path = SEP.join([source_path,genus])
		clusters = os.listdir(genus_path)
		for cluster_ID in clusters:
			seqs_nchars = nchars_jMT_bionj = nchars_mrAIC = nchars_jMT_ml = None
			seqs_ntaxa = ntaxa_jMT_bionj = ntaxa_mrAIC = ntaxa_jMT_ml = None

			logger.info("extracting files for " + genus + " cluster #" + cluster_ID)

			cluster_obj = PhyloObject(genus, cluster_ID)
			if genus + "_" + cluster_ID in source_dict.keys():
				cluster_obj.set_source(*source_dict[genus + "_" + cluster_ID])
			else:
				cluster_obj.set_source("concat", "concat")
			# #parse mrAIC file if exists
			# mraic_dir_location = SEP.join([source_path, genus])
			# mrAIC_values = parse_mrAIC_files.get_mraic_values(mraic_dir_location, cluster_ID)
			# if mrAIC_values is not None:
			# 	mrAIC_AIC_results, mrAIC_AICc_results, mrAIC_BIC_results, ntaxa_mrAIC, nchars_mrAIC = mrAIC_values
			# 	cluster_obj.set_test_results(MR_AIC, AIC, mrAIC_AIC_results)
			# 	cluster_obj.set_test_results(MR_AIC, AICc, mrAIC_AICc_results)
			# 	cluster_obj.set_test_results(MR_AIC, BIC, mrAIC_BIC_results)
			
			#jModelTest for BIONJ tree
			jmt_dir_location = SEP.join([source_path, genus])
			# jmt_bionj_values = parse_jModelTest2_file.get_jmt_values(jmt_dir_location, cluster_ID,
			#                                                    JMODELTEST_BIONJ_OUTPUT_FILENAME)
			# if jmt_bionj_values is not None:
			# 	if jmt_bionj_values == -1:
			# 		logger.warning("problem with parsing jmt_bionj")
			# 	else:
			# 		jmt_AIC_results, jmt_AICc_results, jmt_BIC_results, ntaxa_jMT_bionj, nchars_jMT_bionj = jmt_bionj_values
			# 		cluster_obj.set_test_results(JMODELTEST_BIONJ, AIC, jmt_AIC_results)
			# 		cluster_obj.set_test_results(JMODELTEST_BIONJ, AICc, jmt_AICc_results)
			# 		cluster_obj.set_test_results(JMODELTEST_BIONJ, BIC, jmt_BIC_results)
				
			#jModelTest for ML tree
			jmt_ml_values = parse_jModelTest2_file.get_jmt_values(jmt_dir_location, cluster_ID,
															   JMODELTEST_ML_OUTPUT_FILENAME)
			if jmt_ml_values is not None:
				if jmt_ml_values == -1:
					logger.warning("problem with parsing jmt_ml")
				else:
					jmt_AIC_results, jmt_AICc_results, jmt_BIC_results, ntaxa_jMT_ml, nchars_jMT_ml = jmt_ml_values
					cluster_obj.set_test_results(JMODELTEST_ML, AIC, jmt_AIC_results)
					cluster_obj.set_test_results(JMODELTEST_ML, AICc, jmt_AICc_results)
					cluster_obj.set_test_results(JMODELTEST_ML, BIC, jmt_BIC_results)

			if jmt_ml_values: #or jmt_bionj_values
				#edit msa file path
				msa_file = SEP.join([source_path, genus, cluster_ID, PHYLIP_FILENAME])
				msa_features = parse_MSA_file.get_seqs_values(msa_file, alignment_type)
				cluster_obj.set_msa_features(msa_features)
				tre_file = re.sub(PHYLIP_FILENAME, BIONJ_TREE_FILENAME, msa_file)
				if os.path.exists(tre_file):
					cluster_obj.base_tree_features = parse_newick_tree.parse_tree_file(tre_file)
				phyml_hky_stats_file = SEP.join([source_path, genus, cluster_ID, PHYML_STATS_HKY])
				phyml_hky_tre_file = SEP.join([source_path, genus, cluster_ID, PHYML_TREE_HKY])
				phyml_gtr_stats_file = SEP.join([source_path, genus, cluster_ID, PHYML_STATS_GTR])
				if os.path.exists(phyml_hky_tre_file) and os.path.exists(phyml_gtr_stats_file):
					cluster_obj.set_phyml_stats_features(parse_phyml_for_bionj.get_base_tree_features(
						phyml_hky_stats_file, phyml_gtr_stats_file))
					cluster_obj.phyml_tree_features = parse_newick_tree.parse_tree_file(phyml_hky_tre_file)

				csv_writer.writerow(cluster_obj.get_features_as_list())
		dataFile.flush()

	dataFile.close()
