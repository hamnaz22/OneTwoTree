from helper_functions import get_avg, get_std
import glob
from Bio import Phylo
from Bio.Phylo import Newick
from io import StringIO
import Bio.Nexus.Trees
Tree = Newick.BaseTree.Tree
from pipeline_defs import *

__author__ = 'Shiran'


def clean_non_terminals(tree, clades_dic):
	leaves_dic = {}
	for clade in clades_dic.keys():
		if Tree.is_terminal(clade):
			leaves_dic[clade] = clades_dic[clade]
	return leaves_dic


def extract_tree_properties(tree):
	"""
	returns tree features for object tree (of type Tree)
	"""
	total_dis = Tree.total_branch_length(tree)-1

	#depth of all clades in a dictionary
	clades_level_depth = Tree.depths(tree, unit_branch_lengths=True)
	clades_dist_depth = Tree.depths(tree, unit_branch_lengths=False)

	#keep in the dictionary only terminals
	leaves_level_depth = clean_non_terminals(tree, clades_level_depth)
	leaves_dist_depth = clean_non_terminals(tree, clades_dist_depth)

	#eliminate the +1 of biopython to unrooted trees
	fixed_level_depths = [x-1 for x in leaves_level_depth.values()]
	fixed_dist_depths = [x-1 for x in leaves_dist_depth.values()]

	#avg_level_depth = get_avg(fixed_level_depths)
	avg_dist_depth = get_avg(fixed_dist_depths)

	#std_level_depth = get_std(fixed_level_depths)
	std_dist_depth = get_std(fixed_dist_depths)

	# min_dist_depth = min(fixed_dist_depths)
	# min_level_depth = min(fixed_level_depths)
	#
	max_dist_depth = max(fixed_dist_depths)
	# max_level_depth = max(fixed_level_depths)

	return NewickTreeFeatures(total_dis, avg_dist_depth, std_dist_depth, max_dist_depth)



def parse_tree_string(tree_string):
	#convert tree_string to an io and open it with Phylo.read
	handle = StringIO(tree_string)
	try:
		tree = Phylo.read(handle, "newick")
	except:
		return None

	return extract_tree_properties(tree)


def parse_mraic_tree_file(cluster_path, test_criteria):
	#open tre file
	mraic_tre_file_path = glob.glob(cluster_path + "/*.phy." + test_criteria + "-*.tre")
	if len(mraic_tre_file_path) == 0:
		return [None]*5

	return parse_tree_file(mraic_tre_file_path[0])





def parse_tree_file(tre_file_path):
	#open tre file
	try:
		tree = Phylo.read(tre_file_path, "newick")
	except:
		return None

	return extract_tree_properties(tree)