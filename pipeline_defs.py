from __future__ import division
from sys import platform as _platform
from helper_functions import *
import re

__author__ = 'Shiran'


MSA_999_FILENAME = "seqs-organism.fasta"
MSA_FILENAME = "seqs-organism-concat.fasta"
PHYLIP_FILENAME = "seqs-organism.phy"
DM_PHYLIP_FILENAME = "distance-matrix.phy"
BIONJ_TREE_FILENAME = "bionj_base_tree.tre"

JMT_JAR_FILE = "jModelTest.jar"

MRAIC_TXT_SUFFIX = ".MrAIC.txt"
MRAIC_TXT_FILE = "seqs-organism.phy.MrAIC.txt"
JMODELTEST_BIONJ_OUTPUT_FILENAME = "JMT_BIONJ.output.txt"
JMODELTEST_ML_OUTPUT_FILENAME = "JMT.output.txt"

BASE_TREE_TYPE_BIONJ = "BIONJ"
BASE_TREE_TYPE_ML = "ML"

FASTA_FORMAT = "fasta"
PHYLIP_FORMAT = "phylip-relaxed"

PHYML_STATS_SUFFIX = "_phyml_stats{0}.txt"
PHYML_TREE_SUFFIX = "_phyml_tree{0}.txt"
PHYML_LK_SUFFIX = "_phyml_lk.txt"
PHYML_STATS_OUTPUT = PHYLIP_FILENAME + PHYML_STATS_SUFFIX
PHYML_TREE_OUTPUT = PHYLIP_FILENAME + PHYML_TREE_SUFFIX
PHYML_LK_OUTPUT = PHYLIP_FILENAME + PHYML_LK_SUFFIX

PHYML_STATS_HKY = "base_HKY_PHYML.txt"
PHYML_TREE_HKY = "base_HKY_PHYML.tre"
PHYML_STATS_GTR = "base_GTR_PHYML.txt"
PHYML_TREE_GTR = "base_GTR_PHYML.tre"

PHYML_LNL_SITES_SELECTED_LK = "T1_phyml_lk.txt"
PHYML_LNL_SITES_ML_LK = "T2_phyml_lk.txt"
PHYML_LNL_SITES_SELECTED_STATS = "T1_phyml_stats.txt"
PHYML_LNL_SITES_ML_STATS = "T2_phyml_stats.txt"
PHYML_LNL_SITES_SELECTED_TREE = "T1_phyml_tree.txt"
PHYML_LNL_SITES_ML_TREE = "T2_phyml_tree.txt"
CONSEL_INPUT_FILENAME = "treesLscores.paup"
CONSEL_OUTPUT_FILENAME = "consel_table.txt"

#app types
MR_AIC = "mrAIC"
JMODELTEST = "jMT"
JMODELTEST_ML = "jMT_ML"
JMODELTEST_BIONJ = "jMT_BioNJ"

#statistic tests
AIC = "AIC"
AICc = "AICc"
BIC = "BIC"

#DBs
TRPIPE_DB_NAME = "trpipe"
AA_DB_NAME = "aa_db"
SELETOME_DB_NAME = "selectome"
PANDIT_DB_NAME = "pandit"
RNA_DB_NAME = "RNA"

MODEL_SUBSET = {"JC69": "JC",
				"JC": "JC",
				"K2P": "K2P",
				"K80": "K2P",
				"GTR": "GTR",
				"F81": "F81",
				"HKY": "HKY",
				"SYM": "SYM"}

#PHYML_MODEL_SUBSET = {"JC": "000000",
#				"K2P": "010010",
#				"GTR": "012345",
#				"F81": "000000",
#				"HKY": "010010",
#				"SYM": "012345"}

#MODEL_PARAMETERS = {"JC": 1, "JCI":2, "JCG":2, "JCIG":3,
#					"K2P": 2, "K2PI": 3, "K2PG": 3, "K2PIG": 4,
#					"F81": 4, "F81I": 5, "F81G": 5, "F81IG": 6,
#					"HKY": 5, "HKYI": 6, "HKYG": 6, "HKYIG": 7,
#					"SYM": 6, "SYMI": 7, "SYMG": 7, "SYMIG": 8,
#					"GTR": 9, "GTRI": 10, "GTRG": 10, "GTRIG": 11}

NUM_SUBSTITUTION = {"JC": 1,
				"K2P": 2,
				"GTR": 3,
				"F81": 1,
				"HKY": 2,
				"SYM": 3}


def AICc_penalty(ntaxa, model_params, nchars):
	branch = float(ntaxa)*2-3
	n = float(nchars)
	k = branch+model_params
	try:
		return 2*k*(k+1)/float(n-k-1)
	except:
		return float('NaN')


def break_model_inv_gamma(model):
	inv = gamma = False
	if re.search("^.{2,4}I.?$", model):
		inv = True
	if re.search("^.{2,4}I?G$", model):
		gamma = True
	matrix = MODEL_SUBSET[re.sub("I?G?$", "", model)]

	return [matrix, inv, gamma]


class TestResults:
	def __init__(self, app, test):
		super().__init__()
		self.app = app
		self.test = test
		self.score = None
		self.selected_model = None
		self.base_model = None
		self.transitions_matrix = 0
		self.inv = False     # True/False
		self.gamma = False   # True/False
		self.F = False
		self.tree = NewickTreeFeatures

	def set_values(self, model, score, tree_properties):
		self.score = score
		self.set_model_inv_gamma(model)
		self.set_tree_results(tree_properties)
		return

	def set_model_inv_gamma(self, model):
		[self.base_model, self.inv, self.gamma] = break_model_inv_gamma(model)
		self.selected_model = self.base_model
		if self.inv:
			self.selected_model += "I"
		if self.gamma:
			self.selected_model += "G"
		if self.base_model in ["F81", "HKY", "GTR"]:
			self.F = True
		self.substitution_count = NUM_SUBSTITUTION[self.base_model]

		return

	def set_tree_results(self, newick_tree_results):
		self.tree = newick_tree_results
		return

	def get_features_as_list(self):
		l = ([self.selected_model, self.base_model, self.score, self.inv, self.gamma, self.F, self.substitution_count] + self.tree.get_features_as_list())
		return l


class NewickTreeFeatures:
	def __init__(self, *tree_properties):
		super().__init__()

		total_dis, avg_dist_depth, std_dist_depth, max_depth = tree_properties
		self.max_depth = max_depth
		self.avg_dist_depth = avg_dist_depth
		self.std_dist_depth = std_dist_depth
		self.total_tree_height = total_dis

		return

	def get_features_as_list(self):
		l = ([self.total_tree_height, self.avg_dist_depth, self.std_dist_depth,
			  self.max_depth])
		return l


##class BaseTreeFeatures:
##	def __init__(self):
##		super().__init__()
##
##		self.parsimony = 0
##		self.tree_size = 0
##		self.gamma = 0
##		self.prop_inv = 0
##		self.Ts_ts_ratio = 0
##		self.freq_AC = 0
##		self.freq_AG = 0
##		self.freq_AT = 0
##		self.freq_CG = 0
##		self.freq_CT = 0
##		self.freq_GT = 0
##		self.freqs_std = 0
##		self.matrix_AtoT = 0
##		self.rate = 0
##
##		return
##
##	def get_features_as_list(self):
##		l = [self.parsimony, self.gamma, self.prop_inv, self.Ts_ts_ratio,
##			 self.freq_AC, self.freq_AG, self.freq_AT, self.freq_CG, self.freq_CT, self.freqs_std, self.rate]
##		return l


##class SequencesFeatures:
##	def __init__(self, ntaxa, nchars, alignment_type):
##		super().__init__()
##
##		self.ntaxa = ntaxa
##		self.nchars = nchars
##		self.alignment_type = alignment_type
##
##		self.klogn = self.calculate_klogn()
##		self.gapless_matrix_size = 0
##		self.aicc_penalties = self.calculate_all_AICc_penalty()
##		self.aicc_penalty_deltas = self.calculate_dAICc_penalty_delta(self.aicc_penalties)
##
##		self.gc_content_avg = 0
##		self.gc_content_std = 0
##
##		self.gap_content_avg = 0
##		self.gap_content_std = 0
##
##		self.transition_rate_avg = 0
##		self.transversion_rate_avg = 0
##
##		self.freq_A = 0
##		self.freq_C = 0
##		self.freq_G = 0
##		self.freq_T = 0
##		self.freq_std = 0
##		self.Echanges_per_site = 0
##
##		self.p_invariant_sites = 0
##		self.entropy = 0
##		self.sop_score = 0
##
##	def calculate_klogn(self):
##		k = float(self.ntaxa)*2-3
##		logn = math.log(self.nchars)
##		return  k*logn
##
##	def calculate_dAICc_penalty_delta(self, penalties):
##		d = {}
##		d["p10p7"] = penalties[10]-penalties[7]
##		d["p11p10"] = penalties[11]-penalties[10]
##		d["p10p6"] = penalties[10]-penalties[6]
##		d["p10p5"] = penalties[10]-penalties[5]
##		d["p10p9"] = penalties[10]-penalties[9]
##		d["p10p6"] = penalties[10]-penalties[6]
##		d["p10p5"] = penalties[10]-penalties[5]
##		d["p9p4"] = penalties[9]-penalties[4]
##		d["p9p5"] = penalties[9]-penalties[5]
##
##		return list(d.values())
##
##	def calculate_all_AICc_penalty(self):
##		return [AICc_penalty(self.ntaxa, i, self.nchars) for i in range(12)]
##
##	def set_p_invariant_sites(self, p_invariant_sites):
##		self.p_invariant_sites = p_invariant_sites
##		return
##
##	def set_entropy(self, entropy):
##		self.entropy = entropy
##		return
##
##	def set_nucleotide_frequency(self, nucleotide, frequency):
##		exec("self.freq_" + nucleotide + " = " + str(frequency))
##		return
##
##	def set_gc_values(self, gc_avg, gc_std):
##		self.gc_content_avg = gc_avg
##		self.gc_content_std = gc_std
##		return
##
##	def set_gap_values(self, gap_avg, gap_std):
##		self.gap_content_avg = gap_avg
##		self.gap_content_std = gap_std
##		return
##
##	def set_substitution_rates(self, transition_rate, transversion_rate):
##		self.transition_rate_avg = transition_rate
##		self.transversion_rate_avg = transversion_rate
##		return
##
##	def get_features_as_list(self):
##		l = ([self.alignment_type, self.ntaxa, self.nchars, self.klogn, self.gapless_matrix_size, self.gc_content_avg, self.gc_content_std,
##			  self.gap_content_avg, self.gap_content_std, self.transition_rate_avg,
##			  self.freq_A, self.freq_C, self.freq_G, self.freq_T, self.freq_std, self.Echanges_per_site, self.p_invariant_sites, self.entropy, self.sop_score] +
##			  self.aicc_penalties + self.aicc_penalty_deltas)
##		return l


######lass PhyloObject:

###	def __init__(self, genus, cluster_id):
###		super().__init__()

###		self.id = genus + "_" + cluster_id
###		self.genus = genus
###		self.cluster_id = cluster_id

###		self.source = None
###		self.type = None
###		self.sequences_features = SequencesFeatures

###		self.base_tree_features = NewickTreeFeatures
###		self.phyml_stats_features = BaseTreeFeatures
###		self.phyml_tree_features = NewickTreeFeatures

###		#self.mraic_aic = TestResults
###		#self.mraic_aicc = TestResults
###		#self.mraic_bic = TestResults

###		self.jmt_bionj_aic = TestResults
###		self.jmt_bionj_aicc = TestResults
###		self.jmt_bionj_bic = TestResults

###		self.jmt_ml_aic = TestResults
###		self.jmt_ml_aicc = TestResults
###		self.jmt_ml_bic = TestResults

###		return

###	def set_source(self, c_source, c_type):
###		self.source = c_source
###		self.type = c_type

###	def get_features_as_list(self):
###		# try:
###		# 	mraic_list = self.mraic_aic.get_features_as_list() + self.mraic_aicc.get_features_as_list() + \
###		# 				 self.mraic_bic.get_features_as_list()
###		# except:
###		# 	mraic_list = [None]*11*3

###		try:
###			jmt_bionj_list = self.jmt_bionj_aic.get_features_as_list() + self.jmt_bionj_aicc.get_features_as_list() + \
###					   self.jmt_bionj_bic.get_features_as_list()
###		except:
###			jmt_bionj_list = [None]*11*3

###		try:
###			jmt_ml_list = self.jmt_ml_aic.get_features_as_list() + self.jmt_ml_aicc.get_features_as_list() + \
###					   self.jmt_ml_bic.get_features_as_list()
###		except:
###			jmt_ml_list = [None]*11*3

###		try:
###			seqs_features_list = self.sequences_features.get_features_as_list()
###		except:
###			seqs_features_list = [None]*(19+19) #19 dAICc features

###		#try:
###		#	base_tree_list = self.base_tree_features.get_features_as_list()
###		#except:
###		#	base_tree_list = [None]*7

###		try:
###			phyml_stats_list = self.phyml_stats_features.get_features_as_list()
###		except:
###			phyml_stats_list = [None]*11

###		try:
###			phyml_tree_list = self.phyml_tree_features.get_features_as_list()
###		except:
###			phyml_tree_list = [None]*4

###		l = [self.id] + jmt_bionj_list + jmt_ml_list + seqs_features_list +\
###			phyml_stats_list + phyml_tree_list + [self.source, self.type]
###		return l

###	def set_test_results(self, app, test, results):
###		if app == MR_AIC:
###			if test == AIC:
###				self.mraic_aic = results
###			elif test == AICc:
###				self.mraic_aicc = results
###			elif test == BIC:
###				self.mraic_bic = results
###		elif app == JMODELTEST_BIONJ:
###			if test == AIC:
###				self.jmt_bionj_aic = results
###			elif test == AICc:
###				self.jmt_bionj_aicc = results
###			elif test == BIC:
###				self.jmt_bionj_bic = results
###		elif app == JMODELTEST_ML:
###			if test == AIC:
###				self.jmt_ml_aic = results
###			elif test == AICc:
###				self.jmt_ml_aicc = results
###			elif test == BIC:
###				self.jmt_ml_bic = results
###		return

###	def set_msa_features(self, msa_results):
###		self.sequences_features = msa_results
###		return

###	def set_base_tree_features(self, base_tree_features):
###		self.base_tree_features = base_tree_features
###		return

###	def set_phyml_stats_features(self, phyml_tree_features):
###		self.phyml_stats_features = phyml_tree_features
###		return

###	def convert_to_dictionary(self):
###		#there's a problem with this function... ignore for now
###		d = self.__dict__
###		d.pop("mraic_aic")
###		d.pop("mraic_aicc")
###		d.pop("mraic_bic")
###		d.pop("jmt_aic")
###		d.pop("jmt_aicc")
###		d.pop("jmt_bic")
###		all_items_dic = dict(list(d.items()) + list(self.mraic_aic.__dict__.items()) +
###							 list(self.mraic_aicc.__dict__.items()) + list(self.mraic_bic.__dict__.items()) +
###							 list(self.jmt_bionj_aic.__dict__.items()) + list(self.jmt_bionj_aicc.__dict__.items()) +
###							 list(self.jmt_bionj_bic.__dict__.items()))
###		return all_items_dic

