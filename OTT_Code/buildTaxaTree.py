import datetime
import time
import logging.handlers
import mysql.connector
import os
import re
import json
import zipfile
import argparse
from ctypes import *
from Bio import Phylo
from Bio import AlignIO
from Bio.Nexus import Nexus
from Bio import Alphabet
import random
import pandas as pd

from ott_objects_defs.PloiDbContext import PloiDbContext
from ott_objects_defs.ploidbCommon import *
from ott_objects_defs.ClusterContext import SeqType
from ott_objects_defs.ploidbUtils import *
from ott_objects_defs.MBconfig_Dictionary import mb_block_dict
from ott_objects_defs.handleMultAccessions import handle_multiple_accessions
from ott_objects_defs.createJobFile import *
from ott_objects_defs.PathHelper import PathHelper
from ott_objects_defs.blastSequencesLocaly import *
from ott_objects_defs.getGenusSequences import *
from ott_objects_defs.ITS_usingDB import *

from ott_scripts.Create_Tree import *  # Tree generation functions: part2_create_mrbayes_blk, create_mb_config_file
from ott_scripts.adjustDirFasta import adjust_direction_fasta_using_mafft



"""
This file contains the major methods of ploidb pipeline
"""

use_cached_cluster_results = True


#--------------------------------------------------------------------------------------------------------
#Perform Ultrametric trees for Examl/Raxml without dating, parameter is method: treepl ot dpp
def perform_ML_UltrametricNoDate(context,ml_tree_method,tree_file_path):

	shutil.copyfile(context.summary_dir + '/Result_Tree_' + context.id + '.tre',context.summary_dir +
					'/Result_Tree_' + context.id + '_NotUltra.tre')
	os.chdir(context.xml_dir)
	if ml_tree_method == 'treepl':
		name1,name2 = return_leafEndsNames(tree_file_path)
		msa_conct_length = return_msa_length(context.concat_seqs_report_filename)
		create_NameDate_config(context,msa_conct_length,name1,name2,"1","1")
		cmd_run_tree = 'mpirun -np 4 treePL %s' %context.tree_xml_namedate_config
	elif ml_tree_method == 'dpp':
		# Use R scipt to remove tree labels. Input params: working_dir, tree, tree_out
		tree_out_nolabel = tree_file_path +'_NoLbl'
		scripts_dir = ott_config['general']['OTT_MAIN'] + 'ott_scripts/'
		r_path = ott_config['diff_soft']['R_PATH']
		removeLabels_Rcmd = (r_path + " CMD BATCH '--args working_dir=" + '"' + context.xml_dir + '"' +
							 " tree=" + '"' + tree_file_path + '"' + " tree_out=" + '"' + tree_out_nolabel + '"' + "' "
							 + scripts_dir + "clear_node_labels.R")
		logger.info("Remove labels from tree: %s" % removeLabels_Rcmd)
		os.system(removeLabels_Rcmd)
		# use dos2unix on both fasta and tree files:
		os.system('dos2unix ' + context.concat_seqs_fasta_filename)
		os.system('dos2unix ' + tree_out_nolabel)
		#run dpp on tree file:
		convert_fasta_to_phylip(context.concat_seqs_fasta_filename, context.concat_seqs_fasta_filename+'_phy')
		cmd_run_tree = "mpirun -np 4 dppdiv-mpi-sse3 -in %s -out %s -tre %s -sf 100 " \
					   "-n 10000" %(context.concat_seqs_fasta_filename+'_phy',context.xml_dir+'/UltrametricTree.tree',
									tree_out_nolabel)

	tree_run_out = context.concat_workdir + '/TreeUltra.out'
	tree_run_err = context.concat_workdir + '/TreeUltra.err'
	retval = exec_external_command_redirect_output(cmd_run_tree, tree_run_out, tree_run_err)
	#os.system(cmd_run_tree)
	if retval != 0:
		with open(context.final_status, 'w') as f_status:
			f_status.write("Failed to create Ultra-metric Tree")
		with open(context.working_dir+'/TreeStatusRemark.txt','w') as f_treeRemark:
			f_treeRemark.write("NON Ultra-metric tree")
		raise Exception("Failed to create Ultra-metric Tree")

	#Choose the final tree from the dpp result:
	if ml_tree_method == 'dpp':
		cmd_choose_one = ott_config['diff_soft']['TREE_ANNO'] + ' -heights mean %s.ant.tre outfile.ant.tre.out.treeannotator.tre' %(context.xml_dir+'/UltrametricTree.tree')
		tree_run_out = context.concat_workdir + '/selectOne.out'
		tree_run_err = context.concat_workdir + '/selectOne.err'
		retval = exec_external_command_redirect_output(cmd_choose_one, tree_run_out, tree_run_err)
		if retval != 0:
			with open(context.final_status, 'w') as f_status:
				f_status.write("dpp: Failed to choose the final tree")
			with open(context.working_dir+'/TreeStatusRemark.txt','w') as f_treeRemark:
				f_treeRemark.write("dpp - Failed to choose th final tree")
			raise Exception("dpp: Failed to choose the final tree")
		else:
			shutil.copyfile(context.xml_dir+'/outfile.ant.tre.out.treeannotator.tre',context.summary_dir + '/Result_Tree_' + context.id + '.tre')

	return
#--------------------------------------------------------------------------------------------------------
#Perform Ultrametric trees for Examl/Raxml without dating, parameter is method: treepl ot dpp
def perform_ML_UltrametricWithDate(context,ml_tree_method,tree_file_path,NameDate_missingName_Flag):

	shutil.copyfile(context.summary_dir + '/Result_Tree_' + context.id + '.tre',context.summary_dir + '/Result_Tree_' + context.id + '_NotUltra.tre')
	os.chdir(context.xml_dir)
	if ml_tree_method == 'treepl':
		if NameDate_missingName_Flag == 'No': #verify the names are in the final msa
			cmd_run_tree = 'mpirun -np 4 treePL %s' %context.tree_xml_namedate_config
		else:
			logger.debug("Species given are not on the final tree (since they were not included in final alignment)")
			with open(context.final_status, 'w') as f_status:
				f_status.write("dpp: Failed to choose the final tree")
			with open(context.working_dir+'/TreeStatusRemark.txt','w') as f_treeRemark:
				f_treeRemark.write("dpp - Failed to choose th final tree")
		return
	elif ml_tree_method == 'dpp':
		# Use R scipt to remove tree labels. Input params: working_dir, tree, tree_out
		tree_out_nolabel = tree_file_path +'_NoLbl'
		scripts_dir = ott_config['general']['OTT_MAIN'] + 'ott_scripts/'
		r_path = ott_config['diff_soft']['R_PATH']
		removeLabels_Rcmd = (r_path + " CMD BATCH '--args working_dir=" + '"' + context.xml_dir + '"' +
							 " tree=" + '"' + tree_file_path + '"' + " tree_out=" + '"' + tree_out_nolabel + '"' + "' "
							 + scripts_dir + "clear_node_labels.R")
		logger.info("Remove labels from tree: %s" % removeLabels_Rcmd)
		os.system(removeLabels_Rcmd)
		# use dos2unix on both fasta and tree files:
		os.system('dos2unix ' + context.concat_seqs_fasta_filename)
		os.system('dos2unix ' + tree_out_nolabel)
		#run dpp on tree file:
		convert_fasta_to_phylip(context.concat_seqs_fasta_filename, context.concat_seqs_fasta_filename+'_phy')
		cmd_run_tree = "mpirun -np 4 dppdiv-mpi-sse3 -in %s -out %s -tre %s -sf 100 -n 10000" \
							  %(context.concat_seqs_fasta_filename+'_phy',context.xml_dir+'/UltrametricTree.tree',tree_out_nolabel)

	tree_run_out = context.concat_workdir + '/TreeUltra.out'
	tree_run_err = context.concat_workdir + '/TreeUltra.err'
	retval = exec_external_command_redirect_output(cmd_run_tree, tree_run_out, tree_run_err)
	if retval != 0:
		with open(context.final_status, 'w') as f_status:
			f_status.write("Failed to create Ultra-metric Tree")
		with open(context.working_dir+'/TreeStatusRemark.txt','w') as f_treeRemark:
			f_treeRemark.write("NON Ultra-metric tree")
		raise Exception("Failed to create Ultra-metric Tree")

	#Choose the final tree from the dpp result:
	if ml_tree_method == 'dpp':
		cmd_choose_one = ott_config['diff_soft']['TREE_ANNO'] + ' -heights mean %s.ant.tre outfile.ant.tre.out.treeannotator.tre' %(context.xml_dir+'/UltrametricTree.tree')
		tree_run_out = context.concat_workdir + '/selectOne.out'
		tree_run_err = context.concat_workdir + '/selectOne.err'
		retval = exec_external_command_redirect_output(cmd_choose_one, tree_run_out, tree_run_err)
		if retval != 0:
			with open(context.final_status, 'w') as f_status:
				f_status.write("dpp: Failed to choose the final tree")
			with open(context.working_dir+'/TreeStatusRemark.txt','w') as f_treeRemark:
				f_treeRemark.write("dpp - Failed to choose the final tree")
			raise Exception("dpp: Failed to choose the final tree")
		else:
			#'outfile.ant.tre.out.treeannotator.tre' convert to tre
			nex_chosen = context.xml_dir +'/outfile.ant.tre.out.treeannotator.tre'
			out_Tree = context.summary_dir +'/Result_Tree_%s.tre' % context.id
			Phylo.convert(nex_chosen, 'nexus', out_Tree, 'newick')
			logger.debug("Copy the selected ultrametric tree to the output final tree file")

	return
#--------------------------------------------------------------------------------------------------------
#creates the paragraph for MB config file for Node dating

	#print CONF "constraint split_node = $split_name1 $split_name2;\n";
	#print CONF "calibrate split_node = uniform($split_min_age,$split_max_age);\n";
	#print CONF "prset topologypr = constraints(split_node);\n";
	#print CONF "prset nodeagepr = calibrated;\n";

def create_NodeDate_paragraph_MB(NumOfNodeDates,Name1_list,Name2_list,MinAge_list,MaxAge_list):

	node_idx=0
	mb_node_paragraph=""
	new_str=""
	while node_idx < NumOfNodeDates:
		Node_lbl =  'split_node_' + str(node_idx)
		Name1 = Name1_list[node_idx].replace(" ","_")
		Name2 = Name2_list[node_idx].replace(" ","_")
		new_str = "constraint %s = %s %s;\ncalibrate %s = uniform(%s,%s);\nprset topologypr = constraints(%s);\nprset nodeagepr = calibrated;\n" \
			%(Node_lbl,Name1,Name2,Node_lbl,str(MinAge_list[node_idx]),str(MaxAge_list[node_idx]),Node_lbl)
		node_idx+=1
		mb_node_paragraph+=new_str

	return mb_node_paragraph
#--------------------------------------------------------------------------------------------------------
def add_constraint_blk(context,config_file):

	logger.debug("Combine config file with the constraint Blk generated using R")
	blk_file = context.working_dir + '/concat/constraint_block.txt'
	filenames = [config_file+'_old', blk_file]
	shutil.copyfile(config_file,config_file+'_old')
	with open(config_file, 'w') as outfile:
		for fname in filenames:
			with open(fname) as infile:
				for line in infile:
					if 'end;' in line:
						continue
					else:
						outfile.write(line)
		outfile.write('end;')
	return
#--------------------------------------------------------------------------------------------------------
#This function converts the aligned msa to Nexus format and creates the web partition file (with clusters names)
def make_nexus_msa_and_web_partition(context):
	logger.debug("Performing msa conversion to Nexus format and create web partition file")
	#Convert files to nexus
	Clusters_nexFiles_list = []
	with open(context.concat_workdir + '/fasta-files-to-concat.txt', 'r') as f_list_clusters_files:
		for line in f_list_clusters_files:
			file_name = line.rstrip()
			AlignIO.convert(file_name, 'fasta', file_name + '_nex', 'nexus', alphabet=Alphabet.generic_dna)
			logger.debug("Converted %s to nexus" % (file_name))
			Clusters_nexFiles_list.append(file_name + '_nex')
	if not Clusters_nexFiles_list:	#no selected clusters:
		status_line = 'No clusters found according to user specifications. Try to modify your parameters.'
		logger.debug(status_line)
		with open(context.final_status, 'a') as f_status:
			f_status.write("Failed - %s \n" % status_line)
		raise Exception("Failed - no clusters found")

	logger.debug("Nexus creation, using the following files list:")
	logger.debug(Clusters_nexFiles_list)
	nexi = [(fname, Nexus.Nexus(fname)) for fname in Clusters_nexFiles_list]
	combined = Nexus.combine(nexi)
	combined.write_nexus_data(filename=open(context.concat_workdir + '/MSA_btCOMBINED.nex', 'w'))
	f_web_partition = open(context.summary_dir + '/web_partition.txt', 'w')
	#replace file names with cluster desc:
	f_msa_nexus = open(context.summary_dir + '/msa_nexus.txt', 'w')
	cluster_partition_dict = {}
	with open(context.concat_workdir + '/MSA_btCOMBINED.nex', 'r') as f_src_nex:
		for line in f_src_nex:
			origin_line = line
			flag_written = 'Yes'
			for file_name in Clusters_nexFiles_list:
				if file_name in line:
					for clust_idx in context.loci_description_dict.keys():
						clust_idx_str = '/' + clust_idx + '/'
						if clust_idx_str in file_name:
							line = line.replace(file_name, context.loci_description_dict[clust_idx])
							current_cluster_num = clust_idx
							if 'charpartition combined' in line:
								flag_written = 'remove'
								continue
							flag_written = 'No'
						#take the data needed for web partition
			if flag_written == 'Yes':
				f_msa_nexus.write(origin_line)
			elif flag_written == 'No':
				f_msa_nexus.write(line)
				web_line = line.replace('charset', 'Cluster_%s' % current_cluster_num)
				r = re.compile("'(.*?)'")
				m = r.search(web_line)
				desc_str = m.group(1)
				web_line = web_line.replace("'" + desc_str + "'", '')
				web_line = web_line.replace(';', ', representative annotation: ' + desc_str)
				f_web_partition.write(web_line)
	return

#--------------------------------------------------------------------------------------------------------------
def Create_Tree(context,f_end_status_name):
	Tree_Method = context.UserFlags_dict['Tree_Method']
	Tree_Type = context.UserFlags_dict['Tree_Type'] #Concat or Per cluster/Locus
	if Tree_Type == 'ConcatTree':
		create_Tree_diffMethods(f_end_status_name, Tree_Method, context)
	else:
		#open dir for Trees in summary directory:
		create_dir_if_not_exists(context.summary_dir+'/LocusTrees')
		logger.debug("Open Clusters file: %s" %(context.summary_dir+'/LocusTrees/NumOfClusterTrees.txt'))
		with open(context.summary_dir+'/LocusTrees/NumOfClusterTrees.txt','w') as f_trees_ndexFile:
			#Create Tree per cluster:
			for cluster_context in context.cluster_contexts:
				current_context = context.get_cluster_by_id(cluster_context.cluster_id)
				ClusterID = str(current_context.index)
				#Check if cluster was chosen for the concat file:
				if current_context in context.cluster_contexts_for_concat_tree:
					#for cluster_id in clusters_list:
					create_LocusTree_diffMethods(ClusterID, Tree_Method, context)
					f_trees_ndexFile.write(ClusterID+',')
		#Create zip of all Locus trees:
		LocusTrees_zipName = context.summary_dir + '/LocusTrees_' + context.id + '.zip'
		logger.debug("Locus Trees zip file: %s" %LocusTrees_zipName)
		locusTree_zip = zipfile.ZipFile(LocusTrees_zipName, mode='w')
		for root, dirs, files in os.walk(context.summary_dir+'/LocusTrees'):
			for file in files:
				locusTree_zip.write(os.path.join(root, file), file)
		locusTree_zip.close()
	return
#--------------------------------------------------------------------------------------------------------------
def calc_merged_clusterData(context, cluster_id, line_origin):
	#"ClusterID,Desc,Type,length,SpeciesCnt,SpeciesIds,IncludedInConcat\n")

	x = csv.reader(line_origin)
	split_line = list(x)

	logger.debug(split_line)
	for item in split_line:
		logger.debug(item)

	species_names = []
	desc_origin = ''
	desc_origin = split_line[2][0]
	type_origin = split_line[4][0]
	median_origin = split_line[6][0]
	species_cnt_origin = split_line[8][0]

	count_spc = 0
	with open(context.working_dir + '/concat/' + cluster_id + '/seqs-organism-concat.fasta', 'r') as cluster_alg:
		for line in cluster_alg:
			if '>' in line:
				count_spc += 1
				name_spc = line.replace('>', '')
				name_spc = name_spc.rstrip()
				species_names.append(name_spc)
	spc_list = count_spc

	data_line = """"%s","%s","%s","%s","%d (%s)","%s","yes"\n""" % (
		cluster_id, desc_origin, type_origin, median_origin, count_spc, species_cnt_origin, species_names)
	logger.debug("Merged data line:")
	logger.debug(data_line)
	#""""%s","%s","%s","%s","%d","%s","%s"\n""" %(ClusterID,desc_str,cluster_context.cluster_type_dict[ClusterID],cluster_len,len(cluster_species_list),taxid_accession_list,IncludedInConcat)

	return data_line

def Write_Selected_Clusters_data_File(context, merged_dict):
	data_for_cluster = {}
	logger.debug("Copy selected data to new clusters file:")
	logger.debug(context.rerun_chosen_clusters)
	f_origin_clusters = open(context.UserFlags_dict['OriginJobDir'] + '/SummaryDir/All_Clusters_data.csv', 'r')
	f_selected_clusters = open(context.summary_clusters_data_file, 'w')
	f_selected_clusters.write("ClusterID,Desc,Type,length,SpeciesCnt,SpeciesIds,IncludedInConcat\n")
	#For Merged clusters we need to calculate New data:
	for line in f_origin_clusters:
		logger.debug(line)
		line_split = line.split(',')
		cluster_id = line_split[0].replace('"', '')
		logger.debug("Check if cluster %s was selected" % cluster_id)
		if cluster_id in merged_dict.values():
			data_merged_line = calc_merged_clusterData(context, cluster_id, line)
			f_selected_clusters.write(data_merged_line)
		elif int(cluster_id in context.rerun_chosen_clusters):  #need to copy this line
			f_selected_clusters.write(line)

	return


def edit_file_toMerge(file_to_merge, species_list):
	New_records_to_Merge = []
	with open(file_to_merge, "rU") as toEdit_f:
		for seq_record in SeqIO.parse(file_to_merge, "fasta"):
			organism_name = return_codedName(getPropertyFromFastaSeqHeader(seq_record.description, "organism"))
			logger.debug(organism_name)
			if organism_name not in species_list:
				seq_record.id = organism_name
				seq_record.description = ""
				New_records_to_Merge.append(seq_record)
	f_out = open(file_to_merge + '_edited', 'w')
	SeqIO.write(New_records_to_Merge, f_out, "fasta")

	return


def merge_clusters(cluster_to_merge_idx, main_cluster_idx, context):
	#take all new species from seqs-no-multiple-accessions.fasta in cluster dir and addfrag to concat/1/seqs-organism-concat.fasta of the main cluster (the one we want to merge to)

	#e.g. File to merge: /bioseq/data/results/oneTwoTree/1489995756/7/seqs-no-multiple-accessions.fasta
	#e.g. File to Merge to: /bioseq/data/results/oneTwoTree/1489995756/concat/1/seqs-organism-concat.fasta

	file_to_merge = context.working_dir + '/' + cluster_to_merge_idx + '/seqs-no-multiple-accessions.fasta'
	file_toMergeTo = context.working_dir + '/concat/' + main_cluster_idx + '/seqs-organism-concat.fasta'

	#Edit the file we want to add -> remove species that are already included in the the cluster we merged into and chnage the headers to be the organism alone:
	species_list_inMain = get_organism_from_fasta(file_toMergeTo)
	logger.debug("Species names in file we merged into:")
	logger.debug(species_list_inMain)
	edit_file_toMerge(file_to_merge, species_list_inMain)
	shutil.copyfile(file_toMergeTo, file_toMergeTo + '_origin')
	#addFrag_cmd = 'mafft --addfragments ' + file_to_merge + ' --multipair ' + file_toMergeTo + ' > ' + file_toMergeTo + '_merged'
	addFrag_cmd = 'mafft --quiet --addfragments ' + file_to_merge + '_edited' + ' --multipair ' + file_toMergeTo + ' > ' + file_toMergeTo + '_merged'
	logger.debug(addFrag_cmd)
	os.system(addFrag_cmd)
	shutil.copyfile(file_toMergeTo + '_merged', file_toMergeTo)

	return


def Rerun_flow(rerunJobID, context):
	#Check for taxa editing: get the original final list of species from the previous run and check for updates:
	original_final_species_f = context.UserFlags_dict['OriginJobDir'] + '/SummaryDir/FinalSpeciesList.txt'
	new_species_input_f = context.working_dir + '/userInput.txt'
	update_species_list, context.rerun_remove_species_list, context.rerun_add_species_list = compare_taxa_files(
		original_final_species_f, new_species_input_f)
	#Check for different alignemnt and fasta filter methods (mafft....guidance and such):
	check_methods_diff(context,context.UserFlags_dict['OriginJobDir'])
	logger.debug('Original vs current mafft software is: %s' %context.align_method_diff)
	#rewite taxa_list file:
	with open(context.working_dir + '/taxa_list.txt', 'w') as taxa_list_f:
		for item in update_species_list:
			taxa_list_f.write(item + '\n')
		taxa_list_f.close()
	logger.debug("User edited species list:")
	logger.debug("Remove the following species: %s" % context.rerun_remove_species_list)
	logger.debug("Add new species: %s" % context.rerun_add_species_list)

	#outgroup selection: get the selected outgorup from the original run
	#---------------------------------------------------------------------
	context.outgroupSelection = rerun_get_selected_outgroup(
		os.path.join(context.UserFlags_dict['OriginJobDir'], "OutgroupSelection.txt"))
	shutil.copyfile(context.UserFlags_dict['OriginJobDir'] + "/OutgroupSelection.txt",
					context.working_dir + '/OutgroupSelection.txt')

	#get the constraint file from the original run
	#---------------------------------------------------------------------
	if os.path.exists(context.UserFlags_dict['OriginJobDir']+'/ConstraintTree_user.txt'):
		shutil.copyfile(context.UserFlags_dict['OriginJobDir']+'/ConstraintTree_user.txt',
						context.working_dir + '/ConstraintTree_user.txt')
		constraint_empty='yes'
		with open(context.working_dir+'/ConstraintTree_user.txt') as constraint_f:
			for line in constraint_f:
				if line.strip():    #line is not empty
					constraint_empty='no'
		if constraint_empty=='no':
			context.UserFlags_dict["Constraint_Flag"] = 'On'    #User entered data for Constraint Tree.
			logger.debug("Constraint is on")
		else:
			context.UserFlags_dict["Constraint_Flag"] = 'Off'

	#Get selected clusters for rerun:
	#----------------------------------------------------------------------
	context.Nuc_clusterMode = calc_clustMethod(context.UserFlags_dict['include_Nuc'])
	context.Mt_clusterMode = calc_clustMethod(context.UserFlags_dict['include_Mt'])
	context.Chloro_clusterMode = calc_clustMethod(context.UserFlags_dict['include_CP'])
	logger.debug("Rerun section: Clsuter Type selection")
	logger.debug("Nuc: %s" %context.UserFlags_dict['include_Nuc'])
	logger.debug("Mt: %s" %context.UserFlags_dict['include_Mt'])
	logger.debug("Chloro: %s" %context.UserFlags_dict['include_CP'])
	list_to_include=[]
	if context.Nuc_clusterMode == 1: list_to_include.append('ClustType_NUC')
	if context.Mt_clusterMode == 1: list_to_include.append('ClustType_mtDNA')
	if context.Chloro_clusterMode == 1: list_to_include.append('ClustType_cpDNA')
	logger.debug(context.Nuc_clusterMode)
	logger.debug(context.Mt_clusterMode)
	logger.debug(context.Chloro_clusterMode)
	logger.debug("list_to_include:")
	logger.debug(list_to_include)
	for key in context.UserFlags_dict.keys():
		if 'cluster_' in key:
			if context.UserFlags_dict[key] == 'on':
				#Check Cluster type:
				cluster_index = key.replace('cluster_', '')
				for clustType in list_to_include:
					if os.path.exists(context.UserFlags_dict['OriginJobDir'] + '/' + str(cluster_index) + '/' + clustType):
						context.rerun_chosen_clusters.append(cluster_index)
	if not context.rerun_chosen_clusters:
		logger.info("Rerun - List of Clusters is empty, choose all clusters")
		#Write to final status:
		with open(working_dir + '/SummaryDir/FinalStatus.txt', 'w') as sum_Status_f:
			sum_Status_f.write("No clusters were selected for the rerun, try again")
			return
	else:
		logger.info("Rerun - List of Clusters chosen by the user: %s" % context.rerun_chosen_clusters)

	#Check if need to merge clusters:
	#----------------------------------------------------------------------
	#cluster_1:on  	#Merge_1:4
	for cluster_index in context.rerun_chosen_clusters:
		#if cluster_index in clusters_to_merge_dict.keys():
		#	logger.debug("Skip cluster %s, it was merged to cluster %s" %(cluster_index,clusters_to_merge_dict[cluster_index]))
		#	continue
		#copy clusters dirs into new run:
		src_cluster_dir = context.UserFlags_dict['OriginJobDir'] + '/' + str(cluster_index)
		dest_cluster_dir = context.working_dir + '/'
		copy_dir_to_another_dir(src_cluster_dir, dest_cluster_dir)
		#copy concat clusters dirs:
		src_cluster_dir = context.UserFlags_dict['OriginJobDir'] + '/concat/' + str(cluster_index)
		dest_cluster_dir = context.working_dir + '/concat/'
		copy_dir_to_another_dir(src_cluster_dir, dest_cluster_dir)
		#Enable new alignment for the rerun clusters - only if mafft method is diff from original run:
		curr_cluster_dir = dest_cluster_dir + str(cluster_index)
		if context.align_method_diff == 'true':
			if os.path.exists(curr_cluster_dir+ "/seqs-aligned-concat.fasta"):
				try:
					os.remove(curr_cluster_dir+ "/seqs-aligned-concat.fasta")
				except OSError:
					logger.debug("FAILED to erase old alignment file")
			#Check if ITS cluster:
			cluster_dir = context.working_dir + '/' + str(cluster_index)
			if os.path.exists(cluster_dir +'/ITS_CLUSTER'):
				fasta_for_alignment = cluster_dir + "/seqs-no-multiple-accessions.fasta"
			else:
				fasta_for_alignment = curr_cluster_dir + "/seqs-with-out-group-concat.fasta"
			run_alignment(fasta_for_alignment, curr_cluster_dir+ "/seqs-aligned-concat.fasta",
						  context, curr_cluster_dir)
		if context.UserFlags_dict['FilterMSA_Method'] != 'None':
			msa_file = curr_cluster_dir+ "/seqs-aligned-concat.fasta"
			perform_filter_msa(ploidb_context, msa_file)
		#Check if MSA filter method was chosen (default is none)
		if context.UserFlags_dict['FilterMSA_Method'] != 'None':
			MSA_filter_method = context.UserFlags_dict['FilterMSA_Method']
			logger.debug("Rerun option for MSA filter method is: %s" % MSA_filter_method)


	clusters_to_merge_dict = {}
	merged_clusters_list = []
	for cluster_id in context.rerun_chosen_clusters:
		if context.UserFlags_dict['Merge_' + cluster_id] != 'None':
			merged_clusters_list.append(cluster_id)
			clusters_to_merge_dict[cluster_id] = context.UserFlags_dict['Merge_' + cluster_id]
			merge_clusters(cluster_id, clusters_to_merge_dict[cluster_id], context)

	logger.debug("Clusters that will be merged to others:")
	logger.debug(merged_clusters_list)
	#remove all merged clusters from the list of clusters that will be concatenated:
	for merged_clust in merged_clusters_list:
		context.rerun_chosen_clusters.remove(merged_clust)

	logger.debug("Final list of clusters: after merge")
	logger.debug(context.rerun_chosen_clusters)

	context.cluster_contexts_for_concat_tree = context.rerun_chosen_clusters
	logger.debug("clusters for concat: %s" % context.cluster_contexts_for_concat_tree)

	#Update SummaryDir with all clusters data:
	Write_Selected_Clusters_data_File(context, clusters_to_merge_dict)
	logger.info("Saving clusters data in %s" % context.summary_clusters_data_file)

	#process_seq_and_create_concat_file(context, context.id,context.UserFlags_dict['Outgroup_Flag'])
	#RERUN_process_seq_and_create_concat_file(context, context.id,context.UserFlags_dict['Outgroup_Flag'])

	#Perform updated alignment:
	concat_alignment_rerun(context)
	logger.info("Concat File found: %s" % context.concat_seqs_fasta_filename)
	make_nexus_msa_and_web_partition(context)

	#Result indication file:
	f_end_status_name = context.debug_dir + '/AlignEnd_FinishedExecutionWithConcatFile'
	f_end = open(f_end_status_name, 'w')
	f_end.close()
	#Copy alignment file to Summary dir:
	shutil.copyfile(context.concat_seqs_fasta_filename,
					context.summary_dir + '/' + context.id + "-concat-aligned.fasta")

	#Copy constraint file from original run
	if os.path.exists(context.working_dir+'/ConstraintTree_user.txt'):
		constraint_empty='yes'
		with open(context.working_dir+'/ConstraintTree_user.txt') as constraint_f:
			for line in constraint_f:
				if line.strip():    #line is not empty
					constraint_empty='no'
		if constraint_empty=='no':
			context.UserFlags_dict["Constraint_Flag"] = 'On'    #User entered data for Constraint Tree.
			f_params_web.write('Constraint Tree: ' + context.UserFlags_dict["Constraint_Flag"] + '\n')
		else:
			context.UserFlags_dict["Constraint_Flag"] = 'Off'

	#Run Tree:
	#
	if os.path.exists(f_end_status_name):
		logger.info(
			"======================     Alignment DONE - Start Tree Generation	====================================")
		Tree_Method = context.UserFlags_dict['Tree_Method']
		Tree_Type = context.UserFlags_dict['Tree_Type'] #Concat or Per cluster/Locus
		logger.debug(" Rerun tree type is: %s" % Tree_Type)
		if Tree_Type == 'ConcatTree':
			create_Tree_diffMethods(f_end_status_name, Tree_Method, context)
		else:
			create_dir_if_not_exists(context.summary_dir+'/LocusTrees')
			logger.debug("Open Clusters file: %s" %(context.summary_dir+'/LocusTrees/NumOfClusterTrees.txt'))
			with open(context.summary_dir+'/LocusTrees/NumOfClusterTrees.txt','w') as f_trees_ndexFile:
				#Create Tree per cluster:
				for ClusterID in context.cluster_contexts_for_concat_tree:
					#for cluster_id in clusters_list:
					create_LocusTree_diffMethods(ClusterID, Tree_Method, context)
					f_trees_ndexFile.write(ClusterID+',')
			#Create zip of all Locus trees:
			LocusTrees_zipName = context.summary_dir + '/LocusTrees_' + context.id + '.zip'
			logger.debug("Locus Trees zip file: %s" %LocusTrees_zipName)
			locusTree_zip = zipfile.ZipFile(LocusTrees_zipName, mode='w')
			for root, dirs, files in os.walk(context.summary_dir+'/LocusTrees'):
				for file in files:
					locusTree_zip.write(os.path.join(root, file), file)
			locusTree_zip.close()

	#Enter the number of species in the rerun Alignment:
	#Check for species in fasta file:
	species_num = 0
	with open(context.concat_seqs_fasta_filename, 'r') as f_align:
		for line in f_align:
			if '>' in line:
				species_num += 1
	with open(context.summary_file, 'a') as f_sum:
		f_sum.write("Species in final alignment file: %d\n" % species_num)
	return

def Build_raxml_Tree(Alignment_file, context, xml_model):
	create_dir_if_not_exists(context.xml_dir)
	shutil.copyfile(Alignment_file, context.xml_dir + '/Alignment_file.fasta')
	raxml_out_log = context.xml_dir + '/Raxml_log.txt'
	raxml_out_err = context.xml_dir + '/Raxml_err.txt'

	#raxmlHPC -x 12345 -p 12345 -N 5 -m GTRCAT -s /groups/i.......   -w working_dir -q ....ml_partition -f a

	Create_partitionFile(context.concat_seqs_report_filename, context.concat_xml_partition_filename)
	#Check for outgroup selection:
	if context.UserFlags_dict['Outgroup_Flag'] == 'None':
		if context.UserFlags_dict['Use_BootS_RxOn'] == 'Off':
			raxml_cmd = 'raxmlHPC -s %s -m %s -p 12345 -q %s -w %s -n %s' % (
				Alignment_file, xml_model, context.concat_xml_partition_filename, context.xml_dir, 'raxml_' + context.id)
		else:
			bootstrap_val = str(context.UserFlags_dict['Raxml_BootS_vals'])
			raxml_cmd = 'raxmlHPC -s %s -x 12345 -p 12345 -N %s -m %s -w %s -q %s -n %s -f a' % (
				Alignment_file, bootstrap_val, xml_model, context.xml_dir, context.concat_xml_partition_filename,
				'raxml_' + context.id)
	else:
		genus_outgroup = context.outgroupSelection
		if context.UserFlags_dict['Use_BootS_RxOn'] == 'Off':
			raxml_cmd = 'raxmlHPC -s %s -m %s -o %s -p 12345 -q %s -w %s -n %s' % (
				Alignment_file, xml_model, genus_outgroup, context.concat_xml_partition_filename, context.xml_dir,
				'raxml_' + context.id)
		#raxml_cmd = 'raxmlHPC -m %s -o %s -p 12345 -q %s -s %s -n %s' % (xml_model,genus_outgroup,context.concat_xml_partition_filename,Alignment_file,'raxml_'+context.id)
		else:
			bootstrap_val = str(context.UserFlags_dict['Raxml_BootS_vals'])
			raxml_cmd = 'raxmlHPC -s %s -x 12345 -p 12345 -N %s -m %s -o %s -w %s -q %s -n %s -f a' % (
				Alignment_file, bootstrap_val, xml_model, genus_outgroup, context.xml_dir,
				context.concat_xml_partition_filename, 'raxml_' + context.id)

	if context.UserFlags_dict["Constraint_Flag"] == 'On':
		raxml_cmd += (' -g %s' %(context.working_dir+'/ConstraintTree_user.txt'))


	logger.debug("Execute Raxml: %s" % raxml_cmd)
	os.chdir(os.path.abspath(context.xml_dir))
	exec_external_command_redirect_output(raxml_cmd, raxml_out_log, raxml_out_err)
	#Copy result tree to summary dir:
	if context.UserFlags_dict['Use_BootS_RxOn'] == 'Off':
		shutil.copyfile('RAxML_result.raxml_' + context.id, context.summary_dir + '/Result_Tree_' + context.id + '.tre')
	else:
		shutil.copyfile('RAxML_bipartitions.raxml_' + context.id,
						context.summary_dir + '/Result_Tree_' + context.id + '.tre')

	return


## ----------------------------------------- RAXML Bootstrap Tree
def Build_raxml_Boots_Tree(Alignment_file, context, xml_model,bootstrap_val,outgorup_flag):

	create_dir_if_not_exists(context.xml_dir)
	shutil.copyfile(Alignment_file, context.xml_dir + '/Alignment_file.fasta')
	raxml_out_log = context.xml_dir + '/Raxml_log.txt'
	raxml_out_err = context.xml_dir + '/Raxml_err.txt'

	Create_partitionFile(context.concat_seqs_report_filename, context.concat_xml_partition_filename)
	#Check for outgroup selection:
	if context.UserFlags_dict['Outgroup_Flag'] == 'None':
		raxml_cmd = 'raxmlHPC -s %s -x 12345 -p 12345 -N %s -m %s -w %s -q %s -n %s -f a' % (
			Alignment_file, bootstrap_val, xml_model, context.xml_dir, context.concat_xml_partition_filename,
			'raxml_' + context.id)
	else:
		genus_outgroup = context.outgroupSelection
		raxml_cmd = 'raxmlHPC -s %s -x 12345 -p 12345 -N %s -m %s -o %s -w %s -q %s -n %s -f a' % (
			Alignment_file, bootstrap_val, xml_model, genus_outgroup, context.xml_dir,
			context.concat_xml_partition_filename, 'raxml_' + context.id)

	if context.UserFlags_dict["Constraint_Flag"] == 'On':
		raxml_cmd += (' -g %s' %(context.working_dir+'/ConstraintTree_user.txt'))

	logger.debug("Execute Raxml: %s" % raxml_cmd)
	os.chdir(context.xml_dir)
	exec_external_command_redirect_output(raxml_cmd, raxml_out_log, raxml_out_err)
	#Copy result tree to summary dir:
	shutil.copyfile('RAxML_bipartitions.raxml_' + context.id,
						context.summary_dir + '/Result_Tree_' + context.id + '.tre')
	#shutil.copyfile(context.concat_xml_partition_filename, context.summary_dir + '/partition_file.txt')
	return


# ------------------------------------- RAXML Per Locus TREE ----------------------------------------------
def Build_raxml_LocusTree(Alignment_file, context, xml_model, cluster_id):

	Locus_dir = context.xml_dir + '/' + str(cluster_id) +'/'
	logger.debug("Perform raxml for cluster %s at dir %s" %(cluster_id,Locus_dir))
	create_dir_if_not_exists(Locus_dir)
	shutil.copyfile(Alignment_file, Locus_dir + '/Alignment_file.fasta')
	raxml_out_log = Locus_dir + '/Raxml_log.txt'
	raxml_out_err = Locus_dir + '/Raxml_err.txt'

	Create_partitionFile(context.concat_seqs_report_filename, context.concat_xml_partition_filename)
	#Check for outgroup selection:
	#First Check if outgroup is included in cluster, if not , run tree without outgroup:
	if context.UserFlags_dict['Outgroup_Flag'] != 'None':
		genus_outgroup = context.outgroupSelection
		outgroup_included_flag = is_outgroup_included(genus_outgroup,Alignment_file)
	if context.UserFlags_dict['Outgroup_Flag'] == 'None' or outgroup_included_flag == 'No':
		if context.UserFlags_dict['Use_BootS_RxOn'] == 'Off':
			raxml_cmd = 'raxmlHPC -s %s -m %s -p 12345 -w %s -n %s' % (
				Alignment_file, xml_model, Locus_dir, 'raxml_' + context.id)
		else:
			bootstrap_val = str(context.UserFlags_dict['Raxml_BootS_vals'])
			raxml_cmd = 'raxmlHPC -s %s -x 12345 -p 12345 -N %s -m %s -w %s -n %s -f a' % (
				Alignment_file, bootstrap_val, xml_model, Locus_dir, 'raxml_' + context.id)
	else:
		genus_outgroup = context.outgroupSelection
		if context.UserFlags_dict['Use_BootS_RxOn'] == 'Off':
			raxml_cmd = 'raxmlHPC -s %s -m %s -o %s -p 12345 -w %s -n %s' % (
				Alignment_file, xml_model, genus_outgroup, Locus_dir, 'raxml_' + context.id)
		#raxml_cmd = 'raxmlHPC -m %s -o %s -p 12345 -q %s -s %s -n %s' % (xml_model,genus_outgroup,context.concat_xml_partition_filename,Alignment_file,'raxml_'+context.id)
		else:
			bootstrap_val = str(context.UserFlags_dict['Raxml_BootS_vals'])
			raxml_cmd = 'raxmlHPC -s %s -x 12345 -p 12345 -N %s -m %s -o %s -w %s -n %s -f a' % (
				Alignment_file, bootstrap_val, xml_model, genus_outgroup, Locus_dir,'raxml_' + context.id)

	logger.debug("Execute Raxml: %s" % raxml_cmd)
	os.chdir(Locus_dir)
	exec_external_command_redirect_output(raxml_cmd, raxml_out_log, raxml_out_err)
	#Copy result tree to summary dir:
	if context.UserFlags_dict['Use_BootS_RxOn'] == 'Off':
		shutil.copyfile('RAxML_result.raxml_' + context.id, context.summary_dir+'/LocusTrees/Result_Tree_' + context.id + '_Cluster%s.tre' %(cluster_id))
	else:
		shutil.copyfile('RAxML_bipartitions.raxml_' + context.id, context.summary_dir+'/LocusTrees/Result_Tree_' + context.id + '_Cluster%s.tre' %(cluster_id))


	return


##------------------------------------------------------------------------------------------------------------------
##                          EXAML EXAML EXAML EXAML EXAML EXAML EXAML
##------------------------------------------------------------------------------------------------------------------
def Build_examl_Tree(Alignment_file, context, xml_model):
	#1. run raxml (the regular version - see attached file) which generate (relatively quickly) the following 2 files:
	#      a)  "RAxML_info."  open it and get the number x which is found in the sentence: "Alignment has x distinct alignment patterns"
	#      b)  "RAxML_parsimonyTree."
	#2. parse your alignment to get a binary encoding (examl works with this encoding) with the parse of examl which gets an alignment of phylip foramt, a partitions file and an output name, see attached file.
	#3. run examl (attached) with the binary alignment (the output of the previous step) and the parsimony tree from the first step (you can use another starting tree if you have).
	#to get the number of threads you should divide x by 500 (rule of thumb). you can start without threads if you want.
	#the -S parameter saves memory it is not mandatory.
	#-m parameter is the model examl uses

	create_dir_if_not_exists(context.xml_dir)
	shutil.copyfile(Alignment_file, context.xml_dir + '/Alignment_file.fasta')
	raxml_out_log = context.xml_dir + '/Raxml_log.txt'
	raxml_out_err = context.xml_dir + '/Raxml_err.txt'

	logger.debug("Tree Method: Examl using Raxml for prerun:")
	# 1. raxmlHPC -m GTRGAMMAI -p 12345 -s your_alignment -n output_name -# 1
	os.chdir(context.xml_dir)
	#check if to run Bootstrap
	if context.UserFlags_dict['Outgroup_Flag'] == 'None':
		raxml_cmd = 'raxmlHPC -y -m GTRCAT -p 12345 -s %s -n %s' % (Alignment_file, 'RaxmlPreExaml_' + context.id)  #raxmlHPC -m GTRCAT -p 12345 -s your_alignment -n output_name -# 1
	else:
		genus_outgroup = context.outgroupSelection
		raxml_cmd = 'raxmlHPC -y -m GTRCAT -p 12345 -s %s -o %s -n %s' % (Alignment_file, genus_outgroup, 'RaxmlPreExaml_' + context.id)  #raxmlHPC -m GTRCAT -p 12345 -s your_alignment -n output_name -# 1

	#Regular flow
	logger.debug("Execute prerun for Examl: %s" % raxml_cmd)
	exec_external_command_redirect_output(raxml_cmd, raxml_out_log, raxml_out_err)

	info_fileName = context.xml_dir + '/RAxML_info.RaxmlPreExaml_' + context.id
	f_info = open(info_fileName, 'r')
	for line in f_info:
		if 'distinct alignment patterns' in line:
			Xnum = re.split('Alignment has | distinct alignment patterns', line)[1]
			if int(Xnum) / 500 < 1:
				Xthred = 1
			else:
				Xthred = min(int(int(Xnum) / 500),8)
			logger.debug(
				"Raxml prerun: X number is %s, Xthred is %d ('Alignment has | distinct alignment patterns')" % (
					Xnum, Xthred))
	# 2. Convert fasta to phylip format
	# 2.1 Create partition file
	output_phy_file = Alignment_file + '_phy'
	convert_fasta_to_phylip(Alignment_file, output_phy_file)
	Create_partitionFile(context.concat_seqs_report_filename,
						 context.concat_xml_partition_filename)  #The input is the list od clusters in the alignment file

	# 3. parse-examl -s your_alignment_in_phylip_format -q partitions_file -m DNA -n choose_a_name
	parse_out_name = context.id + '_parseExaml_out'
	os.chdir(context.xml_dir)
	parse_cmd = 'parse-examl -s %s -q %s -m DNA -n %s' % (
		output_phy_file, context.concat_xml_partition_filename, parse_out_name)
	logger.debug("Execute parse cmd for Examl: %s" % parse_cmd)
	os.system(parse_cmd)
	#exec_external_command_redirect_output(parse_cmd,raxml_out_log,raxml_out_err)

	# 4. module load rocks-openmpi
	# 5. mpirun -np num_threads examl -s binary_alignment -t parsimony_tree -m GTRGAMMAI -n output_name -S
	examl_cmd = 'mpirun -np %d examl -s %s -t %s -m %s -n %s' % (
		Xthred, parse_out_name + '.binary', 'RAxML_parsimonyTree.' + 'RaxmlPreExaml_' + context.id, xml_model,
		'Examl_' + context.id)

	#Constraint tree options:
	#First we need to check if the constarint tree includes all species in the fional alignment, otherwise it will not work:
	if context.UserFlags_dict["Constraint_Flag"] == 'On':
		if check_for_constriant_species(context) == 'Yes':
			examl_cmd += (' -g %s' %(context.working_dir+'/ConstraintTree_user.txt'))
		else:
			with open(context.working_dir+'/TreeStatusRemark.txt','w') as f_treeRemark:
				f_treeRemark.write("The constraint tree was not applied since some of its taxa were missing from the alignment")

	logger.debug("Execute Examl cmd: %s" % examl_cmd)
	#Create .sh file for module load + final run:
	sh_file_examl = open(context.xml_dir + '/examl_run.sh', 'w')
	sh_file_examl.write('#!/bin/tcsh\n')
	sh_file_examl.write('\n')
	sh_file_examl.write('cd %s/\n' % context.xml_dir)
	sh_file_examl.write('module load rocks-openmpi\n')
	sh_file_examl.write('%s\n' % examl_cmd)
	sh_file_examl.close()
	os.system('tcsh %s' % (context.xml_dir + '/examl_run.sh'))

	shutil.copyfile('ExaML_result.Examl_' + context.id,
					context.summary_dir + '/Result_Tree_' + context.id + '.tre')
	return

##----------------------------------------
## EXAML Bootstrap !!!
##----------------------------------------
def Build_examl_Boots_Tree(Alignment_file, context, xml_model,bootstrap_val,outgorup_flag):

	create_dir_if_not_exists(context.xml_dir)
	shutil.copyfile(Alignment_file, context.xml_dir + '/Alignment_file.fasta')
	raxml_out_log = context.xml_dir + '/Raxml_log.txt'
	raxml_out_err = context.xml_dir + '/Raxml_err.txt'

	xml_align_file = 'Alignment_file.fasta'
	#xml_align_file = 'Alignment_file.fasta.reduced'
	xml_partition_file = context.id+'_xml_partition'

	logger.debug("Tree Method: Examl (using Raxml to create bootstrap parsimony trees):")
	os.chdir(context.xml_dir)

	# 1. Create input phylip.binary
	Alignment_file_phy = 'Alignment_file.fasta_phy'
	logger.debug("Create phylip alignment file: %s" %(Alignment_file_phy))
	convert_fasta_to_phylip(xml_align_file, Alignment_file_phy)
	Create_partitionFile(context.concat_seqs_report_filename,xml_partition_file)

	## 2. parse-examl -s your_alignment_in_phylip_format -q partitions_file -m DNA -n choose_a_name
	parse_name = 'phy_bin'
	parse_cmd = 'parse-examl -s %s -q %s -m DNA -n %s' % (Alignment_file_phy, xml_partition_file, parse_name)
	logger.debug("Execute parse cmd to transform the alignment file into a binary file: %s" % parse_cmd)
	os.system(parse_cmd)

	# 3. raxmlHPC -N 14 -b 12345 -f j -m GTRGAMMA -s Alignment_file.fasta -n _BooRep (output will be INPUT.BSxx (e.g. Alignment_file.fasta.BS0)
	bootstrapRep_name =  '_BooRep'
	reduced = '.reduced'
	#Check if reduced files exists -> means that the alignment file contains undefined positions:

	if context.UserFlags_dict['Outgroup_Flag'] == 'None':
		raxml_reduced_cmd = 'raxmlHPC -y -m GTRCAT -p 12345 -s %s -q %s -n pars' %(xml_align_file,xml_partition_file)
	else:
		genus_outgroup = context.outgroupSelection
		raxml_reduced_cmd = 'raxmlHPC -y -m GTRCAT -p 12345 -s %s -q %s -o %s -n pars' %(xml_align_file,xml_partition_file,genus_outgroup)


	os.system(raxml_reduced_cmd)
	if os.path.exists(xml_align_file+reduced):
		xml_align_file+=reduced
	if os.path.exists(xml_partition_file+reduced):
		xml_partition_file+=reduced
	boots_rep_cmd = 'raxmlHPC-AVX -N %s -b 12345 -f j -m GTRGAMMA -s %s -q %s -n %s' % (bootstrap_val,xml_align_file,xml_partition_file,bootstrapRep_name)
	logger.debug("Execute bootsrap using raxml: %s" % boots_rep_cmd)
	os.system(boots_rep_cmd)

	info_fileName = context.xml_dir + '/RAxML_info._BooRep'
	f_info = open(info_fileName, 'r')
	for line in f_info:
		if 'distinct alignment patterns' in line:
			Xnum = re.split('Alignment has | distinct alignment patterns', line)[1]
			if int(Xnum) / 500 < 1:
				Xthred = 1
			else:
				Xthred = int(int(Xnum) / 500)
			logger.debug("Raxml prerun: X number is %s, Xthred is %d ('Alignment has | distinct alignment patterns')" % (Xnum, Xthred))

	# 4. Compute starting trees for each BS :
	# raxmlHPC-AVX -y -s ../2_BSrep/$phy.phylip.BS$r -m GTRCAT -n vc_$i\_BS$r -p 12$i$r
	files_to_concat=[]
	for Boots_idx in range(0, int(bootstrap_val)):
		os.chdir(context.xml_dir)
		BS_idx_Dir = (context.xml_dir+'/BS_'+str(Boots_idx))
		partition_BS = xml_partition_file +'.BS'+str(Boots_idx)
		alignment_BS = xml_align_file+'.BS'+str(Boots_idx)
		create_dir_if_not_exists(BS_idx_Dir)
		shutil.move(alignment_BS, BS_idx_Dir+'/')
		shutil.move(partition_BS, BS_idx_Dir+'/')
		os.chdir(BS_idx_Dir)
		random_number = random.randint(1, 1000000)

		parsimony_cmd = 'raxmlHPC-AVX -y -s %s -q %s -m GTRCAT -p %d -n _ParsTree%d ' %(alignment_BS,partition_BS,random_number,Boots_idx)
		logger.debug("Execute parsimony Tree cmd for each Bootstrap: %s" % parsimony_cmd)
		os.system(parsimony_cmd)

	# Perform parse to binary for each tree: parse-examl -s Alignment_file.fasta.reduced.BS9 -q 1498381169_xml_partition.reduced.BS9 -m DNA -n ParseBinar
	# Michal: this part was added since Gblocks runs failed...probably since they caused more seqs to be similar. then the initial binary file can't be used since
		parse_eachTree_cms = 'parse-examl -s %s -q %s -m DNA -n ParseBinar' %(alignment_BS,partition_BS)
		logger.debug("Execute parse-examl for the reduced files: %s" % parse_eachTree_cms)
		os.system(parse_eachTree_cms)

	#5. Run ExaML with mpi
		parse_BST_parsTree = 'RAxML_parsimonyTree._ParsTree'+str(Boots_idx)
		#examl_bs_tree_cmd = 'mpirun -np %d examl -s ../phy_bin.binary -n examlOut -m GAMMA -t %s' %(Xthred,parse_BST_parsTree)
		examl_bs_tree_cmd = 'mpirun -np %d examl -s ParseBinar.binary -n examlOut -m GAMMA -t %s' %(Xthred,parse_BST_parsTree)

		#Constraint tree options:
		#First we need to check if the constarint tree includes all species in the fional alignment, otherwise it will not work:
		if context.UserFlags_dict["Constraint_Flag"] == 'On':
			if check_for_constriant_species(context) == 'Yes':
				examl_cmd += (' -g %s' %(context.working_dir+'/ConstraintTree_user.txt'))


		logger.debug("Execute examl run per BS parsimony input tree: %s" % examl_bs_tree_cmd)
		os.system(examl_bs_tree_cmd)
		#Create .sh file for module load + final run:

		##sh_file_examl = open(BS_idx_Dir + '/examl_run_BS'+str(Boots_idx)+'.sh', 'w')
		##sh_file_examl.write('#!/bin/tcsh\n')
		##sh_file_examl.write('\n')
		##sh_file_examl.write('#$ -p -1\n')
		##sh_file_examl.write('#$ -l h=!(compute-7-1|compute-8-13)\n')
		##sh_file_examl.write('\n')
		##sh_file_examl.write('cd %s/\n' % BS_idx_Dir)
		##sh_file_examl.write('module load rocks-openmpi\n')
		##sh_file_examl.write('%s\n' % examl_bs_tree_cmd)
		##sh_file_examl.close()
		##os.system('qsub examl_run_BS'+str(Boots_idx)+'.sh') #examl_run_BS0.sh


		###for filename in glob.glob(BS_idx_Dir+'/*binaryCheckpoint*'):
		#	os.remove(filename)
		files_to_concat.append(BS_idx_Dir+'/ExaML_result.examlOut')

	#Need to perform poling for a few hours to make sure all were done and only then concat into file:
	if polling_BS_Done(context.xml_dir,bootstrap_val) == 0:
		with open(context.xml_dir+'/TREESall', 'w') as outfile:
			for fname in files_to_concat:
				with open(fname) as infile:
					outfile.write(infile.read())
	else:
		logger.debug("Polling on BS dirs failed !!!!")
		return
		#concat_cmd = 'cat %s  > %s' %(parse_BST_parsTree,context.xml_dir+'/'+TREES_file)
	#6. Concatenate all trees :
	#Move all Raxml previous info files into a different directory
	create_dir_if_not_exists(context.xml_dir+'/ExamlInfo/')
	os.system('mv RAxML_info*.* %s' %(context.xml_dir+'/ExamlInfo/'))
	#Then run the final cmd for Bootstrap tree:
	os.chdir(context.xml_dir)
	#exml_bootstrap_tree_cmd = 'raxmlHPC -f b -m %s -s phy_bin.binary -z TREESall -t %s -n Bootstrap_tree' % (xml_model,context.xml_dir+'/BS_0/ExaML_result.examlOut')
	exml_bootstrap_tree_cmd = 'raxmlHPC -f b -m GTRGAMMA -s phy_bin.binary -z TREESall -t %s -n Bootstrap_tree' % (context.xml_dir+'/BS_0/ExaML_result.examlOut')
	if (outgorup_flag != 'None'): exml_bootstrap_tree_cmd += ' -o %s' %context.outgroupSelection
	logger.debug("Execute raxml for Examl Final Bootstrap Tree: %s" % exml_bootstrap_tree_cmd)
	os.system(exml_bootstrap_tree_cmd)

	shutil.copyfile('RAxML_bipartitions.Bootstrap_tree',context.summary_dir + '/Result_Tree_' + context.id + '.tre')
	return


def Build_PhyML_Tree(f_end_status_name, context):
	# Covert fasta file to phylip:
	concat_dir = context.working_dir + "/concat/"
	concat_fasta = concat_dir + context.id + "-concat-aligned.fasta"
	concat_fasta_phylip = concat_dir + context.id + "-concat-aligned_phy.phy"
	if os.path.exists(concat_fasta):
		convert_fasta_to_phylip(concat_fasta, concat_fasta_phylip)
		phyml_cmd = "PhyML -i " + concat_fasta_phylip
		logger.info("Running tree cmd: %s " % phyml_cmd)
		os.system(phyml_cmd)
		return
	else:
		logger.info(
			"=======================================================================================================")
		logger.info(
			"===========      Tree Generation Failed - missing fasta concat file !!!!!     =========================")
		logger.info(
			"=======================================================================================================")
		return


def create_trees_allFormats(context, outgrouup_selection, concat_dir):
	scripts_dir = ott_config['general']['OTT_MAIN'] + 'ott_scripts/'
	r_path = ott_config['diff_soft']['R_PATH']
	#Create map tree and parse trees:
	AllTreesR_command = (
		r_path+ " CMD BATCH '--args working_dir=" + '"' + concat_dir + '"' + " scripts_dir=" + '"' + scripts_dir +
		'"' + " id=" + '"' + context.id + '"' + " OUT_GROUP=" + '"' + outgrouup_selection + '"' + "' " +
		scripts_dir + "run_mb_pipe.R " + concat_dir +
		"run_mb_pipe.Rout " + "r" + context.id)
	logger.info(" running: %s" % AllTreesR_command)
	os.system(AllTreesR_command)
	map_newick = concat_dir + "/parsemb_map_tree.tre"
	trees_newick = concat_dir + "/parsemb_trees.tre"
	con_nexus = concat_dir + "/mb.out.con.tre"
	#Copy alignment file to Summary dir:
	logger.debug("Copy tree files to Summary dir and create nexus and newicl format")
	shutil.copyfile(map_newick, context.summary_dir + '/parsemb_map_tree.tre')
	shutil.copyfile(trees_newick, context.summary_dir + '/parsemb_trees.tre')
	shutil.copyfile(con_nexus, context.summary_dir + '/mb.out.con.tre')
	shutil.copyfile(concat_dir + "/mb_config.nex", context.summary_dir + '/mb_config.nex')
	#New
	map_nexus = context.summary_dir + "/parsemb_map_tree.nex.tre"
	con_newick = context.summary_dir + "/mb.out.con.newick.tre"
	Phylo.convert(map_newick, 'newick', map_nexus, 'nexus')
	Phylo.convert(con_nexus, 'nexus', con_newick, 'newick')
	logger.debug("Newick/Nexus files are ready at: %s" % context.summary_dir)
	#Create Result_Tree.tre:
	#f_con_tree = open(con_newick,'r')
	f_map_tree = open(map_newick, 'r')
	for line in f_map_tree:
		with open(context.summary_dir + '/Result_Tree_' + context.id + '.tre', 'w') as f_tree_summary:
			f_tree_summary.write(line)
			f_tree_summary.close()
	logger.info("========================================= THE END ==========================================")
	#remove all mb out files, jmt job files:
	mb_list_to_remove = glob.glob(concat_dir + '/mb.out*')
	for file in mb_list_to_remove:
		os.remove(file)
	jmt_job_files = glob.glob(concat_dir + '/*jmt-*')
	for file in jmt_job_files:
		os.remove(file)
	return


def Build_MB_Tree(f_end_status_name, context, nodeDateTxt):
	#This means all data for alignment is ready -> copy files to Summary dir:
	#species_names_file.csv

	outgrouup_selection = context.UserFlags_dict['Outgroup_Flag']
	#Continue to Tree generation if AlignEnd_FinishedExecutionWithConcatFile exists:
	concat_dir = context.working_dir + "/concat/"
	genus = context.id
	mb_outputfilePath = concat_dir
	#Check if tree file exists:
	if os.path.exists(concat_dir + '/mb.out.con.tre'):
		create_trees_allFormats(context, outgrouup_selection, concat_dir)
		return
	if (f_end_status_name == context.debug_dir + '/AlignEnd_FinishedExecutionWithConcatFile'):
		fastaFilesList = (concat_dir + "fasta-files-to-concat.txt")
		if os.path.exists(fastaFilesList):
			if context.UserFlag_mbUserModel != 'None':
				mb_UserModel_part12(fastaFilesList,context.UserFlag_mbUserModel)
			else:
				JMTjobs_dirs_list = part1_create_jmt(context.working_dir, context.ott_scripts_dir, fastaFilesList, context.standAlone_Flag)
		else:
			logger.info(
				"=======================================================================================================")
			logger.info(
				"===========  Tree Generation Failed - missing file fasta-files-to-concat.txt  =========================")
			logger.info(
				"=======================================================================================================")
			return
		if context.UserFlag_mbUserModel == 'None':
			#Check periodicly if Jmodel was completed for all files and then run MBayes:
			if (call_Check_Jmt_Done(JMTjobs_dirs_list) != 0):
				f_jmt = open(context.working_dir + "/JModel_Failed", 'w')
				f_jmt.close()
				logger.info(
					"=======================================================================================================")
				logger.info(
					"===========  Polling Time for JmodelTest has expired, Check jobs status....   =========================")
				logger.info(
					"=======================================================================================================")
				statusFail_LogFile(context, 'Polling Time for JmodelTest has expired')
				return
		TreeFile = (context.working_dir + "/concat/mb.out.con.tre")
		if os.path.exists(TreeFile):
			logger.info("===========    Tree file exists at " + context.working_dir + " !!!   ===========")
		else:
			#Create MB config file for MBayes run:
			fastaFilesList = (concat_dir + "fasta-files-to-concat.txt")
			#Create MB Block According to Model:
			if context.UserFlag_mbUserModel == 'None':  #Model will be set according to user selection and not JmodelTest:
				if (part2_create_mrbayes_blk(genus, fastaFilesList)) == 1:
					logger.debug("    ***   FAILED get_Nparams_ForModel !!!   *** ")
					return
			#Node Dating params - in case we have the paragraph for the config file:
			if nodeDateTxt != 'Empty':
				mb_params_list = [context.UserFlags_dict['ngen'], context.UserFlags_dict['nchains'],
								  context.UserFlags_dict['samplefreq'], context.UserFlags_dict['burninFrac'],
								  context.UserFlags_dict['checkFreq'],
								  nodeDateTxt]
				logger.debug("Found Node Dating txt to add to config file:")
				logger.debug(nodeDateTxt)
				split_flag='yes'
			else:
				split_flag='no'
				mb_params_list = [context.UserFlags_dict['ngen'], context.UserFlags_dict['nchains'],
				  context.UserFlags_dict['samplefreq'], context.UserFlags_dict['burninFrac'],
				  context.UserFlags_dict['checkFreq']]

			#Make sure clock vairables are passed correctly:
			Relaxed_options_list = ['cpp', 'tk02', 'igr']
			relaxed_branch_option_param = 'None'
			#Check if User chose 'UltraMetric Tree' option, then set the clock model to be Uniform:
			if "UltraMetricFlag" in context.UserFlags_dict:
				if context.UserFlags_dict["UltraMetricFlag"] == "Yes":
					clock_model_param = 'uniform'
				else:
					clock_model_param = context.UserFlags_dict['clock_model']
			else:
				clock_model_param = context.UserFlags_dict['clock_model']
			if clock_model_param in Relaxed_options_list:
				relaxed_branch_option_param = context.UserFlags_dict['relaxed_branch_option']
			logger.debug("split_flag")
			logger.debug(split_flag)
			create_mb_config_file(genus, context.working_dir, fastaFilesList, mb_params_list,
								  context.ott_scripts_dir, clock_model_param, relaxed_branch_option_param,split_flag)

		mb_config_File = (concat_dir + "mb_config.nex")
		if os.path.exists(mb_config_File):
			#Check if Constraint Tree was given so the blk will be added to the config file:
			if context.UserFlags_dict["Constraint_Flag"] == 'On':
				add_constraint_blk(context,mb_config_File)
			# Run MrBayes job:
			#check if tree exists:
			if not os.path.exists(concat_dir + 'mb.out.con.tre'):
				mb_command = "mb " + mb_outputfilePath + "mb_config.nex" + " > " + mb_outputfilePath + "MB_LOG.txt "
			else:
				mb_command = ''
			if outgrouup_selection != 'None':
				genus_outgroup = context.outgroupSelection
				if genus_outgroup == 'NULL':
					logger.info(
						"=======================================================================================================")
					logger.info(
						"===========  FAILED - No Outgroup selection !!!						       =========================")
					logger.info(
						"=======================================================================================================")
					statusFail_LogFile(context, 'No Outgroup selection !!!')
					return
				if genus_outgroup == 'FailedToRead_Outgroup':
					logger.info(
						"=======================================================================================================")
					logger.info(
						"===========  FAILED to read outgroup file !!!							       =========================")
					logger.info(
						"=======================================================================================================")
					statusFail_LogFile(context, 'FAILED to read outgroup file !!!')
					return
			#job_filename = createJobFile.create_job_file(genus + "_MrBayes\n", mb_command, genus + "_MrBayes.sh", concat_dir)
			mb_config_jobidFile = open(context.working_dir + "/concat/mb_JobID.txt", 'w')
			logger.debug(mb_command)
			os.system(mb_command)
			#exec_external_command_redirect_output(mb_command,context.working_dir + "/concat/MB_OE.txt",context.working_dir + "/concat/MB_ERR.txt")
			mb_config_jobidFile.close()
			#Check if tre file exists:
			if os.path.exists(concat_dir + 'mb.out.con.tre'):
				shutil.copyfile(concat_dir + '/MB_LOG.txt', context.summary_dir + '/MB_LOG.txt')
				create_trees_allFormats(context, outgrouup_selection, concat_dir)
				logger.info(
					"=======================================================================================================")
				logger.info(
					"===========                            Tree Created Successfully !!!                        ===========")
				logger.info(
					"=======================================================================================================")
			#with open(context.final_status,'a') as f_status:
			#	f_status.write("PASS - Tree Created Successfully\n")
			else:
				logger.info(
					"=======================================================================================================")
				logger.info(
					"===========            Tree criation FAILED, Check output file for detailes     !!!         ===========")
				logger.info(
					"=======================================================================================================")
				with open(context.final_status, 'w') as f_status:
					f_status.write("MrBayes Tree generation FAILED, check input files and parameters")
				shutil.copyfile(concat_dir + '/MB_LOG.txt', context.summary_dir + '/MB_LOG.txt')
				return
		else:
			logger.info(
				"===========  \n\n" + genus + " No mb_config.nex file !!! check for errors   =========================")
			with open(context.final_status, 'w') as f_status:
				f_status.write("MrBayes Failed - configuration file is missing")
				return 'fail'
	else:
		logger.info(
			"=======================================================================================================")
		logger.info("===========  Alignment Failed - %s!!!  =========================" % f_end_status_name)
		logger.info(
			"=======================================================================================================")
		with open(context.final_status, 'a') as f_status:
			f_status.write("Alignment Failed - No tree generation\n")
		return


#This function will take the clusters chosen by the user and create new alignment files:
def concat_alignment_rerun(context):
	logger.info("Starting Rerun concatenation")
	create_dir_if_not_exists(context.mrbayes_concat_work_dir)

	indexFileName_Origin = context.UserFlags_dict['OriginJobDir'] + '/concat/fasta-files-to-concat.txt'
	indexFileName = context.fasta_files_list_to_concat
	f_origin_concat_list = open(indexFileName_Origin, 'r')
	with open(indexFileName, "w") as handle:
		for line in f_origin_concat_list:
			r = re.compile('concat/(.*?)/')
			logger.info("Searching for cluster number r -> %s" % r)
			m = r.search(line)
			logger.info("Searching for cluster number m -> %s" % m)
			if m:
				clust_num = m.group(1)
				if clust_num in context.rerun_chosen_clusters:
					logger.info("cluster %s is included in rerun" % clust_num)
					new_line = line.replace(context.UserFlags_dict['OriginJobDir'], context.working_dir)
					logger.debug(line)
					logger.debug(new_line)
					handle.write(new_line)

	#run_dos2unix(indexFileName)
	concatAlignPAth = ott_config['general']['OTT_MAIN'] + '/ott_scripts'
	concat_align_command = "perl " + concatAlignPAth + "/ConcateAlignments.pl %s %s %s NO NA NA" % \
						   (indexFileName, context.concat_seqs_fasta_filename, context.concat_seqs_report_filename)
	exec_external_command_redirect_output(concat_align_command, context.concat_log_out_filename,
										  context.concat_log_err_filename)


#This is the version used before concat_alignment (Michal D.) !!!
def concat_alignment_SingleOutgroup(context):
	dict_clusterMode = {0: "Don't Include Clusters", 1: "Use NonConcat Clusters"} #, 2: 'Use Concat Clusters'}
	logger.info("Starting concatenation")
	create_dir_if_not_exists(context.mrbayes_concat_work_dir)

	NucClusterMode = context.Nuc_clusterMode
	MtClusterMode = context.Mt_clusterMode
	ChloroClusterMode = context.Chloro_clusterMode

	logger.debug("Clusters Mode: Nuc,Mt,Chloroplast")
	logger.debug("Nuc: " + dict_clusterMode[NucClusterMode])
	logger.debug("Mt: " + dict_clusterMode[MtClusterMode])
	logger.debug("Chloroplast: " + dict_clusterMode[ChloroClusterMode])

	indexFileName = context.fasta_files_list_to_concat
	with open(indexFileName, "w") as handle:
		#This section will decide which clusters to concat: default is to use all clusters separately:
		if NucClusterMode == 1:
			for cluster_context in context.get_concat_clusters_by_type(SeqType.Nuc):  #MD changed SeqType.Nuc
				aligned_fasta_filname = cluster_context.all_seqs_for_concat_tree_creation_fasta_filename
				logger.debug("Writing %s to %s" % (aligned_fasta_filname, indexFileName))
				handle.write(aligned_fasta_filname + "\n")
		elif (context.Nuc_cluster is not None) and (NucClusterMode == 2):
			logger.debug("Adding %s to concat " % context.concat_by_type_fastas[SeqType.Nuc])
			handle.write(context.concat_by_type_fastas[SeqType.Nuc] + "\n")

		if MtClusterMode == 1:
			for cluster_context in context.get_concat_clusters_by_type(SeqType.Mt):  #MD changed SeqType.Nuc
				aligned_fasta_filname = cluster_context.all_seqs_for_concat_tree_creation_fasta_filename
				logger.debug("Writing %s to %s" % (aligned_fasta_filname, indexFileName))
				handle.write(aligned_fasta_filname + "\n")
		elif (context.Mt_cluster is not None) and (MtClusterMode == 2):
			logger.debug("Adding %s to concat " % context.concat_by_type_fastas[SeqType.Mt])
			handle.write(context.concat_by_type_fastas[SeqType.Mt] + "\n")

		if ChloroClusterMode == 1:
			for cluster_context in context.get_concat_clusters_by_type(SeqType.Chloroplast):  #MD changed SeqType.Nuc
				aligned_fasta_filname = cluster_context.all_seqs_for_concat_tree_creation_fasta_filename
				logger.debug("Writing %s to %s" % (aligned_fasta_filname, indexFileName))
				handle.write(aligned_fasta_filname + "\n")
		elif (context.chloroplast_cluster is not None) and (ChloroClusterMode == 2):
			logger.debug("Adding %s to concat " % context.concat_by_type_fastas[SeqType.Chloroplast])
			handle.write(context.concat_by_type_fastas[SeqType.Chloroplast] + "\n")


	#run_dos2unix(indexFileName)
	concatAlignPAth = ott_config['general']['OTT_MAIN'] + '/ott_scripts'
	concat_align_command = "perl " + concatAlignPAth + "/ConcateAlignments.pl %s %s %s NO NA NA" % \
						   (indexFileName, context.concat_seqs_fasta_filename, context.concat_seqs_report_filename)
	exec_external_command_redirect_output(concat_align_command, context.concat_log_out_filename,
										  context.concat_log_err_filename)

#-----------------------------------------------------------------------------------------------------------------
def Create_ClusterType_dict(context):
	cluster_idx_list = context.rerun_chosen_clusters
	cluster_id_Type_dict = {}
	for idx in cluster_idx_list:
		cluster_dir = context.working_dir + '/' + idx + '/'
		if os.path.exists(cluster_dir + 'ClustType_cpDNA'):
			cluster_id_Type_dict[idx] = 'cpDNA'
		elif os.path.exists(cluster_dir + 'ClustType_NUC'):
			cluster_id_Type_dict[idx] = 'NUC'
		elif os.path.exists(cluster_dir + 'ClustType_mtDNA'):
			cluster_id_Type_dict[idx] = 'mtDNA'

	logger.debug("Rerun Cluster Id vs Type:")
	logger.debug(cluster_id_Type_dict)
	return cluster_id_Type_dict

def is_seqs_for_concat_tree_exists(context):
	for cluster_context in context.cluster_contexts_for_concat_tree:
		if not os.path.exists(cluster_context.all_seqs_for_concat_tree_creation_fasta_filename):
			logger.debug("File not exists %s " % cluster_context.all_seqs_for_concat_tree_creation_fasta_filename)
			return False

	return True


def calc_clustMethod(int_1_str):
	if int_1_str == 'off':
		return 0
	elif int_1_str == 'on':
		return 1

def Kill_jmt_jobs(context):
	#Kill all JMT jobs
	#qstat | grep -w OTT-CarEpi

	jmt_jobs_list_f = context.mrbayes_concat_work_dir + '/jmt_JobsList.txt'
	out_f = open(context.mrbayes_concat_work_dir + '/jmt_Status_atEND.txt', 'a')
	with open(jmt_jobs_list_f, 'r') as f_jobs:
		for line in f_jobs:
			jmt_job_name = line.strip()
			os.system('qstat | grep -w %s > %s' % (jmt_job_name, out_f))
	return


def out_function_OTT(context):
	#This function should clean all jobs/files that should be removed in case
	# OTT failed or stopped (no matter what was the reason)
	logger.debug("Performing Cleanup, exit OTT !!!")
	Kill_jmt_jobs(context)
	return

# Update all User Flags dict context variables from params file:
def UpdateContextFlags(context):
	f_flags = open(context.working_dir + '/params.txt', 'r')
	for line in f_flags:
		if 'date:' in line:
			split_line = line.split('date:')
			context.UserFlags_dict['date'] = split_line[1]
			logger.debug("Flag: date is: %s" % split_line[1])
			continue
		if line[0] == '#' or (':' not in line):
			continue
		line = line.strip()
		split_line = line.split(':')
		flag_name = split_line[0]
		flag_val = split_line[1].strip()  #split(' ')[0]
		context.UserFlags_dict[flag_name] = flag_val
		logger.debug("Flag: %s is: %s" % (flag_name, flag_val))

	context.Nuc_clusterMode = calc_clustMethod(context.UserFlags_dict['include_Nuc'])
	context.Mt_clusterMode = calc_clustMethod(context.UserFlags_dict['include_Mt'])
	context.Chloro_clusterMode = calc_clustMethod(context.UserFlags_dict['include_CP'])

	#In case of a rerun: copy all flags for species and outgroup from the original params file:
	if context.UserFlags_dict['OriginJobID'] != 'None':
		logger.debug("Load outgroup params from original job dir")
		f_origin_flags = open(context.UserFlags_dict['OriginJobDir'] + '/params.txt', 'r')
		for line in f_origin_flags:
			if line[0] == '#' or (':' not in line):
				continue
			line = line.strip()
			split_line = line.split(':')
			flag_name = split_line[0]
			flag_val = split_line[1].strip()  #split(' ')[0]
			if flag_name == 'Outgroup_Flag':
				logger.debug("%s: %s" % (flag_name, flag_val))
				context.UserFlags_dict['Outgroup_Flag'] = flag_val
			elif flag_name == 'Outgroup_User':
				logger.debug("%s: %s" % (flag_name, flag_val))
				context.UserFlags_dict['Outgroup_User'] = flag_val
				#Check if user entered species name,in case its a high ranked taxa send failure notice:
				if ' ' not in flag_val:
					status_line = 'User outgroup name must be a species and not high ranked taxa'
					logger.debug(status_line)
					with open(context.final_status, 'a') as f_status:
						f_status.write("Failed - %s \n" % status_line)
					raise Exception("Failed - User outgroup name must be a species and not high ranked taxa")

	logger.info(context.UserFlags_dict)

	#prepare file to present on OTT web:
	f_params_web = open(context.working_dir + '/params_for_web.txt', 'w')
	f_params_web.write('Submitted on: ' + context.UserFlags_dict['date'] + '\n')

	f_params_web.write('Species selection options:\n')
	f_params_web.write('Include Species Descendants: ' + context.UserFlags_dict['SpeciesDescendants'] + '\n')
	f_params_web.write('Name resolution: ' + context.UserFlags_dict['NameResType'] + '\n')
	f_params_web.write('Database selection: ' + context.UserFlags_dict['NR_DB_name'] + '\n')

	f_params_web.write('Filter options:\n')
	f_params_web.write('Filter Hybrids: ' + context.UserFlags_dict['Filter_Hybrids'] + '\n')
	f_params_web.write('Filter SubSp: ' + context.UserFlags_dict['Filter_SubSp'] + '\n')
	f_params_web.write('Filter Unresolved: ' + context.UserFlags_dict['Filter_Unresolved'] + '\n')
	f_params_web.write('Merge SubSp/Variants: ' + context.UserFlags_dict['Merge_Subsp'] + '\n')

	if context.UserFlags_dict['OriginJobID'] == 'None':
		f_params_web.write('\nRerun is Off\n')
	else:
		f_params_web.write('\nRerun is On\n')
		f_params_web.write('OriginJobID: ' + context.UserFlags_dict['OriginJobID'] + '\n')
	f_params_web.write('\nOutgroup_Flag: ' + context.UserFlags_dict['Outgroup_Flag'] + '\n')
	if 'Outgroup_User' in context.UserFlags_dict:
		if context.UserFlags_dict['Outgroup_Flag'] == 'User':
			#Check if outgroup user is legal:
			user_out = context.UserFlags_dict['Outgroup_User'].strip()
			if ' ' not in user_out:
				status_line = 'User outgroup name must be a species and not high ranked taxa'
				logger.debug(status_line)
				with open(context.final_status, 'a') as f_status:
					f_status.write("Failed - %s \n" % status_line)
				raise Exception("Failed - User outgroup name must be a species and not high ranked taxa")
			else:
				f_params_web.write('Outgroup_User: ' + user_out + '\n')

	if context.UserFlags_dict['include_Nuc'] != 'on' or context.UserFlags_dict['include_Mt'] != 'on' \
			or context.UserFlags_dict['include_CP'] != 'on':
		f_params_web.write('Genome selection:\n')
		f_params_web.write('include_Nuc: ' + context.UserFlags_dict['include_Nuc'] + '\n')
		f_params_web.write('include_Mt: ' + context.UserFlags_dict['include_Mt'] + '\n')
		f_params_web.write('include_CP: ' + context.UserFlags_dict['include_CP'] + '\n')
	else:
		f_params_web.write('Genome selection: all\n')

	f_params_web.write('\nClustering options:\n')
	f_params_web.write('Clustering method: ' + context.UserFlags_dict['ClusteringMethod'] + '\n')
	if context.UserFlags_dict['ClusteringMethod'] == 'Ortho':
		f_params_web.write('Seq_Ratio: ' + context.UserFlags_dict['orthoSeq_Ratio'] + '\n')
		f_params_web.write('Ortho_Inflation: ' + context.UserFlags_dict['Ortho_Inflation'] + '\n')
	if context.UserFlags_dict['ClusteringMethod'] == 'BlastClust':
		f_params_web.write('Precent Identity: ' + context.UserFlags_dict['BC_percentIdentity'] + '\n')
		f_params_web.write('Coverage length: ' + context.UserFlags_dict['BC_CoverageLength'] + '\n')

	f_params_web.write('\nAlignment parameters:\n')
	f_params_web.write('MSA_Software: ' + context.UserFlags_dict['MSA_Software'] + '\n')
	if context.UserFlags_dict['MSA_Software'] == 'MAFFT':
		f_params_web.write('MAFFT_maxiterate: ' + context.UserFlags_dict['MAFFT_maxiterate'] + '\n')
		f_params_web.write('PairwiseAlignmentMethod: ' + context.UserFlags_dict['PairwiseAlignmentMethod'] + '\n')
		f_params_web.write('\nFilterMSA_Method: ' + context.UserFlags_dict['FilterMSA_Method'] + '\n')
	if context.UserFlags_dict['FilterMSA_Method'] == 'Trimal':
		f_params_web.write('Trimal_CutOff: ' + context.UserFlags_dict['Trimal_CutOff'] + '\n')
	elif context.UserFlags_dict['FilterMSA_Method'] == 'GUIDANCE':
		f_params_web.write('Guidance_RowCol: ' + context.UserFlags_dict['Guidance_RowCol'] + '\n') #Both, RowsOnly, ColOnly

	f_params_web.write('\nPhylogeny inference tool:\n')
	f_params_web.write('Tree_Method: ' + context.UserFlags_dict['Tree_Method'] + '\n')
	f_params_web.write('Tree_Type: ' + context.UserFlags_dict['Tree_Method'] + '\n')
	if context.UserFlags_dict['Tree_Method'] == 'MrBayes':
		f_params_web.write('User_MrBayes_Model: ' + context.UserFlags_dict['User_MrBayes_Model'] + '\n')
		f_params_web.write('clock_model: ' + context.UserFlags_dict['clock_model'] + '\n')
		f_params_web.write('ngen: ' + context.UserFlags_dict['ngen'] + '\n')
		f_params_web.write('relaxed_branch_option: ' + context.UserFlags_dict['relaxed_branch_option'] + '\n')
		f_params_web.write('samplefreq: ' + context.UserFlags_dict['samplefreq'] + '\n')
		f_params_web.write('nchains: ' + context.UserFlags_dict['nchains'] + '\n')
		f_params_web.write('burninFrac: ' + context.UserFlags_dict['burninFrac'] + '\n')
		f_params_web.write('checkFreq: ' + context.UserFlags_dict['checkFreq'] + '\n')
	if context.UserFlags_dict['Tree_Method'] == 'RAxML':
		f_params_web.write('User_RAxML_Model: ' + context.UserFlags_dict['User_RAxML_Model'] + '\n')
		f_params_web.write('Use_BootS_RxOn: ' + context.UserFlags_dict['Use_BootS_RxOn'] + '\n')
		if context.UserFlags_dict['Use_BootS_RxOn'] != 'off':
			f_params_web.write('Raxml_BootS_vals: ' + context.UserFlags_dict['Raxml_BootS_vals'] + '\n')
	if context.UserFlags_dict['Tree_Method'] == 'ExaML':
		f_params_web.write('User_ExaML_Model: ' + context.UserFlags_dict['User_ExaML_Model'] + '\n')
		f_params_web.write('Use_BootS_ExOn: ' + context.UserFlags_dict['Use_BootS_ExOn'] + '\n')
		if context.UserFlags_dict['Use_BootS_ExOn'] != 'off':
			f_params_web.write('Examl_BootS_vals: ' + context.UserFlags_dict['Examl_BootS_vals'] + '\n')

	#Check for Constraint file/data:
	#In case Not OTT site - create an empty constraint file:
	if os.path.exists(context.working_dir+'/ConstraintTree_user.txt'):
		constraint_empty='yes'
		with open(context.working_dir+'/ConstraintTree_user.txt') as constraint_f:
			for line in constraint_f:
				if line.strip():    #line is not empty
					constraint_empty='no'
		if constraint_empty=='no':
			context.UserFlags_dict["Constraint_Flag"] = 'On'    #User entered data for Constraint Tree.
			f_params_web.write('Constraint Tree: ' + context.UserFlags_dict["Constraint_Flag"] + '\n')
		else:
			context.UserFlags_dict["Constraint_Flag"] = 'Off'
	else:
		#Create an empty file for constained tree:
		f_const = open(context.working_dir+'/ConstraintTree_user.txt','w')
		context.UserFlags_dict["Constraint_Flag"] = 'Off'
		f_const.close

	return

def initLogger(log_filename, debug_log_filename, error_log_filename):
	logger.setLevel(logging.DEBUG)
	logger.handlers = []

	# Main logger - for normal usage
	normal_handler = logging.FileHandler(log_filename)
	formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
	normal_handler.setFormatter(formatter)
	normal_handler.setLevel(logging.INFO)
	logger.addHandler(normal_handler)
	# Roll over on application start
	# logger.handlers[0].doRollover()

	# For debugging/tracing
	detailedHandler = logging.FileHandler(debug_log_filename)
	detailedHandler.setFormatter(formatter)
	detailedHandler.setLevel(logging.DEBUG)
	logger.addHandler(detailedHandler)
	# logger.handlers[1].doRollover()

	# For Errors/warnings
	errorHandler = logging.FileHandler(error_log_filename)
	errorHandler.setFormatter(formatter)
	errorHandler.setLevel(logging.WARNING)
	logger.addHandler(errorHandler)


#
# Call orthomcl in order to cluster the family sequences. These are the parameters for the script
# (1) fasta file contains all the sequences,
# (2) tables names suffix (unique),
# (3) percent match cutoff (for blast) - default: 50,
# (4) inflation parameter (for mcl) - default: 1.5
# (5) directory contains all the scripts
# (6) GenBank file contains the sequences of this genus
#
def cluster_sequences_file_with_orthomcl(fasta_seq_filename, id, scripts_path, gb_seq_filename, output_dir,
										 context):
	if (context.its_support):
		logger.info("ITS is ON - running ITS clustering")

		# For the non ITS clusters
		fasta_without_its = os.path.join(output_dir, context.id + "-allseq-no-its.fasta")
		remove_its_command = "perl " + scripts_path + "getSeqsWithoutITS.pl " + fasta_seq_filename + " " + fasta_without_its
		logger.info("Calling getSeqsWithoutITS (for ITS removal) - " + remove_its_command)

		exec_external_command_redirect_output(command_to_exec=remove_its_command,
											  outfile=output_dir + "/getSeqsWithoutITS.out",
											  errfile=output_dir + "/getSeqsWithoutITS.err")

		fasta_file_to_split = fasta_without_its
		# This is for the ITS clusters
		cluster_ITS(context, fasta_seq_filename, context.cluster_script_dir + "/clustering/", gb_seq_filename)

	else:
		logger.info("ITS is OFF - will not cluster use specific ITS clustering")
		fasta_file_to_split = fasta_seq_filename

	# TODO: currently split by feature is only working for non ITS sequences. Is this the correct behavior?

	if (context.split_by_feature):
		fasta_splitted_by_feature = fasta_seq_filename + "-splitted-by-feature.fasta"
		split_by_feature_command = "perl " + scripts_path + "splitByFeature_GB.pl " + fasta_file_to_split + " " \
								   + fasta_splitted_by_feature + " " + gb_seq_filename
		logger.info("Calling splitByFeature_GB (splitting the sequences by feature ) - " + split_by_feature_command)
		exec_external_command_redirect_output(command_to_exec=split_by_feature_command,
											  outfile=output_dir + "/splitByFeature.out",
											  errfile=output_dir + "/splitByFeature.err")
		orthomcl_fasta_input = fasta_splitted_by_feature
	else:
		orthomcl_fasta_input = fasta_file_to_split

	#Get Orthomcl parameters:
	orthomcl_inflation = context.UserFlags_dict['Ortho_Inflation']

	#check if input for orthomcl has enough data: (michal: need to check how to handle small input to orthomcl)
	#number_of_taxa,list_of_taxa=count_taxa_in_fasta(orthomcl_fasta_input)
	#if (number_of_taxa >= 5):
	if context.UserFlags_dict['ClusteringMethod'] == 'Ortho':
		filter_ratio = context.UserFlags_dict['orthoSeq_Ratio']
		orthomcl_bin_path = ott_config['orthomcl']['orthomcl_bin']
		logger.debug("Clustering parameters: Method is Orthomcl, sequence ratio = %s, Inflation Index = %s" % (filter_ratio, orthomcl_inflation))
		clusteringCommand = "perl " + scripts_path + "orthomcl.pl " + orthomcl_fasta_input + " " + id + " 50 " \
						  + orthomcl_inflation + " " + scripts_path + " " + gb_seq_filename + " " + output_dir + " " \
						  + args.config_filename + " " + filter_ratio + " " + orthomcl_bin_path
		logger.info("Calling orthomcl - " + clusteringCommand)
	else: #perform BlastClust
		PrecentIdentity = context.UserFlags_dict['BC_percentIdentity']
		CoverageLength = context.UserFlags_dict['BC_CoverageLength']
		logger.debug("Clustering paramters: Methods is BlastClust, PrecentIdentity = %s, CoverageLength = %s, Inflation Index = %s" % (PrecentIdentity, CoverageLength, orthomcl_inflation))
		convert_BlastClust_output_toGroups(context,PrecentIdentity,CoverageLength,context.working_dir+'/clustering/'+id+'-allseq-no-its.fasta')
		clusteringCommand = "perl " + scripts_path + "blast_clust.pl " + orthomcl_fasta_input + " " + id + " 50 " \
						  + orthomcl_inflation + " " + scripts_path + " " + gb_seq_filename + " " + output_dir + " " \
						  + args.config_filename + " " + PrecentIdentity  #" 50 1.5 " + scripts_path + " " + gb_seq_filename + " " + output_dir + " " + args.config_filename
		logger.info("Calling BlastClust - " + clusteringCommand)
	#execute the selected Clustering method:
	exec_external_command_redirect_output(command_to_exec=clusteringCommand, outfile=output_dir + "/orthoMCL.out",
										  errfile=output_dir + "/orthoMCL.err")

	if os.path.exists(context.clustering_results_dir):
		shutil.rmtree(context.clustering_results_dir)
	logger.debug("Copying %s to %s" % (output_dir + "_seqs", context.clustering_results_dir))
	shutil.copytree(output_dir + "_seqs", context.clustering_results_dir)

	logger.info("DONE Clustering sequences data for " + id + " using %s" %context.UserFlags_dict['ClusteringMethod'])


def FilterLongSeqs_ITSFile(in_fasta_file, out_fasta_file):
	speciesSeqListToWrite = list()
	seq_length_dict = {}
	f_out = open(out_fasta_file, 'w')
	Total_length = 0
	Total_seq_num = 0
	for seq_record in SeqIO.parse(in_fasta_file, "fasta"):
		Total_seq_num += 1
		seq_length_dict[seq_record.id] = len(seq_record)
		Total_length += len(seq_record)
	AVG_length_1p5 = 1.5 * ( Total_length / Total_seq_num)
	AVG_length_cuttOff = 2 * ( Total_length / Total_seq_num)
	if sum(seq_length_dict[seq_record.id] < AVG_length_1p5 for seq_record in
		   SeqIO.parse(in_fasta_file, "fasta")) / Total_seq_num > 0.9:
		AVG_length_cuttOff = AVG_length_1p5
	for seq_record in SeqIO.parse(in_fasta_file, "fasta"):
		if seq_length_dict[seq_record.id] < AVG_length_cuttOff:
			speciesSeqListToWrite.append(seq_record)
	SeqIO.write(speciesSeqListToWrite, f_out, "fasta")
	print(Total_seq_num)
	#f_out.write("%s: %f\n" % (seq_record.id,seq_length_dict[seq_record.id]))

	return

def cluster_ITS(context, fasta_seq_filename, scripts_path, gb_seq_filename):
	output_dir = context.cluter_its_dir

	# Filtering the non ITS sequences
	fasta_only_its = os.path.join(output_dir, context.id + "-allseq-its-only_.fasta")
	fasta_only_its_filtered = os.path.join(output_dir, context.id + "-allseq-its-only.fasta")

	remove_none_its_command = "perl " + scripts_path + "getSeqsWithITS.pl " + fasta_seq_filename + " " + fasta_only_its
	logger.info("Calling getSeqsWithITS (for ITS clustering) - " + remove_none_its_command)
	retval = exec_external_command_redirect_output(command_to_exec=remove_none_its_command,
												   outfile=output_dir + "/getSeqsWithITS.out",
												   errfile=output_dir + "/getSeqsWithITS.err")

	if retval != 0:
		raise Exception("Failed to get only ITS seqs. Aborting")

	#Check if Large species had ITS sequences and Add them:
	if os.path.exists(context.largeSpeciesITS):
		if os.stat(context.largeSpeciesITS).st_size == 0:
			logger.info("Large species did not include ITS sequences")
		else:
			logger.info("Large species include ITS sequences, adding to main its file")
			its_f = open(fasta_only_its, 'a')
			with open(context.largeSpeciesITS, 'rU') as itsLarge_file:
				seq_records = SeqIO.parse(itsLarge_file, 'fasta')
				SeqIO.write(seq_records, its_f, 'fasta')
			its_f.close()

	#In case there is #
	#if os.stat(fasta_only_its).st_size == 0:
	#	logger.debug("No ITS sequences for this run")
	#	return
	#In case there is #
	if os.stat(fasta_only_its).st_size == 0:
		logger.debug("No ITS sequences for this run")
		return

	#Check for species number limit : min_species_in_cluster
	min_species_in_cluster = int(ott_config['general']['min_species_in_cluster'])
	list_its_taxas=[]
	with open(fasta_only_its,'r') as f_its_only:
		for line in f_its_only:
			if '>' in line:
				exp = r"taxonid\|(.*?)(\||\\n)"
				m = re.search(exp, line)
				if m:
					taxonId = m.group(1)
					logger.debug("ITS seq found for taxId %s" %taxonId)
					if taxonId not in list_its_taxas:
						list_its_taxas.append(taxonId)
	if len(list_its_taxas) < min_species_in_cluster:
		logger.debug("ITS sequences are available for less than %s species, there will be no ITS cluster in this case")
		context.its_min_taxa = 'no'
		return

	scripts_dir = context.cluster_script_dir + "/clustering"
	#MD: add filter for sequences > 2*adg seq length in ITS file:
	FilterLongSeqs_ITSFile(fasta_only_its, fasta_only_its_filtered)

	#Insert an option for ITS flow using the new Database:
	logger.debug("Using ITS python version - DEBUG mode:")
	main_ITS_py(context,output_dir, fasta_only_its_filtered, gb_seq_filename, scripts_dir, args.config_filename, context.should_run_guidance, context.UserFlags_dict['MSA_Software'])


	#Call filter method if needed:
	if context.UserFlags_dict['FilterMSA_Method'] == 'Trimal':
		logger.debug("ITS MSA Filter method: Trimal")
		runTrimal(output_dir + '/oneSeqPerSpecies.msa', context.UserFlags_dict['Trimal_CutOff'])
	elif context.UserFlags_dict['FilterMSA_Method'] == 'Gblocks':
		logger.debug("ITS MSA Filter method: Gblocks")
		runGblocks(output_dir + '/oneSeqPerSpecies.msa', context.UserFlags_dict)


def get_schema_name_for_genus(genus_taxon_id):
	return "o_%s" % genus_taxon_id


def drop_orthomcl_scehma(genus_taxon_id):
	logger.info("Dropping orthomcl schema for genus %s" % genus_taxon_id)

	db_host = ott_config['ott_mysql']['hostname']
	db_user = ott_config['ott_mysql']['username']
	db_password = ott_config['ott_mysql']['password']

	db = mysql.connector.connect(host=db_host, user=db_user, password=db_password)

	logger.debug("Connected to DB " + db_host)
	cur = db.cursor()
	try:
		logger.debug("getting scehma name for %s" % genus_taxon_id)
		schema_name = get_schema_name_for_genus(genus_taxon_id)
		logger.debug("scehma name for %s is %s" % (genus_taxon_id, schema_name))
		logger.debug("Dropping database %s" % schema_name)
		cur.execute("DROP DATABASE %s;" % schema_name)
		logger.debug("database dropped")
	finally:
		cur.close()
		db.close()


# Make sure that orthomcl has enough permissions
# GRANT ALL ON %s.* to %s@'10.1.%%' IDENTIFIED BY 'mysql11';
def create_orthomcl_scehma(genus_taxon_id):
	logger.info("Creating orthomcl schema for genus %s" % genus_taxon_id)

	db_host = ott_config['ott_mysql']['hostname']
	db_user = ott_config['ott_mysql']['username']
	db_password = ott_config['ott_mysql']['password']

	db = mysql.connector.connect(host=db_host, user=db_user, password=db_password)

	logger.debug("Connected to DB " + db_host)
	cur = db.cursor()
	logger.info("Connected to DB")

	try:
		schema_name = get_schema_name_for_genus(genus_taxon_id)
		logger.info("About to create DB %s" % schema_name)
		cur.execute("CREATE DATABASE %s;" % schema_name)
	finally:
		cur.close()
		db.close()


def drop_orthomcl_tables(genus_taxon_id):
	db_host = ott_config['ott_mysql']['hostname']
	db_user = ott_config['ott_mysql']['username']
	db_password = ott_config['ott_mysql']['password']

	db = mysql.connector.connect(host=db_host, user=db_user, password=db_password)
	cur = db.cursor()

	try:
		cur.execute(" DROP TABLE IF EXISTS SimilarSequences" + genus_taxon_id + ";")
		cur.execute(" DROP TABLE IF EXISTS InParalog" + ";")
		cur.execute(" DROP TABLE IF EXISTS Ortholog" + ";")
		cur.execute(" DROP TABLE IF EXISTS CoOrtholog" + ";")
		cur.execute(" DROP TABLE IF EXISTS BestInterTaxonScore" + ";")
		cur.execute(" DROP TABLE IF EXISTS BestQueryTaxonScore" + ";")
		cur.execute(" DROP VIEW IF EXISTS InterTaxonMatch" + ";")
		cur.execute(" DROP TABLE IF EXISTS SimilarSequences" + ";")

		cur.execute(" DROP TABLE IF EXISTS InParalog" + genus_taxon_id + ";")
		cur.execute(" DROP TABLE IF EXISTS Ortholog" + genus_taxon_id + ";")
		cur.execute(" DROP TABLE IF EXISTS CoOrtholog" + genus_taxon_id + ";")
		cur.execute(" DROP TABLE IF EXISTS BestInterTaxonScore" + genus_taxon_id + ";")
		cur.execute(" DROP TABLE IF EXISTS BestQueryTaxonScore" + genus_taxon_id + ";")
		cur.execute(" DROP TABLE IF EXISTS UniqueSimSeqsQueryId" + genus_taxon_id + ";")
		cur.execute(" DROP VIEW IF EXISTS InterTaxonMatch" + genus_taxon_id + ";")
	finally:
		cur.close()
		db.close()


def run_guidance(context, fasta_filename, output_dir, dataset="MSA"):
	logger.info("STARTING Guidance execution on %s. Work dir is %s" % (fasta_filename, output_dir))

	if context.UserFlags_dict['MSA_Software'] == 'ClustalOmega':
		guidance_command = 'perl %s --seqFile %s ' \
						   '--msaProgram CLUSTALO --clustalo %s --seqType nuc --outDir %s --dataset %s --program GUIDANCE --TreeAlg FastTree --Tree_Param \\\\-fastest' % (
							   ott_config['diff_soft']['Guidance'],fasta_filename, ott_config['diff_soft']['ClastaO'], output_dir, dataset)
	else:
		guidance_command = 'perl %s --program GUIDANCE --seqFile %s ' \
						   '--msaProgram MAFFT --seqType nuc --outDir %s --dataset %s --MSA_Param "\-\-adjustdirection"' % (
							ott_config['diff_soft']['Guidance'],fasta_filename, output_dir, dataset)

	exec_external_command_redirect_output(command_to_exec=guidance_command, outfile=output_dir + "/guidance.out",
										  errfile=output_dir + "/guidance.err")

	logger.info("FINISHED Guidance execution on %s" % fasta_filename)


def run_msa_clastal(local_dir, fasta_filename, fasta_final_output):
	MSA_cmd = '%s -i %s -o %s --outfmt=fasta' % (ott_config['diff_soft']['ClastaO'],fasta_filename, fasta_final_output)
	# Save headers for reconstruction
	header_dict = {}
	FullDesc_records = []
	seq_idx = 0
	for seq_record in SeqIO.parse(fasta_filename, "fasta"):
		header_dict[seq_idx] = seq_record.description
		seq_idx += 1

	exec_external_command_redirect_output(command_to_exec=MSA_cmd, outfile=local_dir + "/mafft_log.out",
										  errfile=local_dir + "/mafft_log.err")

	return


def run_alignment(fasta_filename, fasta_final_output, context, local_dir):
	logger.info("Guidance option is Off, running only Alignment on: %s" % fasta_filename)
	#Check Alignment Flag options:
	if context.UserFlags_dict['MSA_Software'] == 'MAFFT':
		MSA_cmd = 'mafft --nuc --quiet --adjustdirection %s > %s' % (fasta_filename, fasta_final_output)
		exec_external_command_redirect_output(command_to_exec=MSA_cmd, outfile=local_dir + "/mafft_log.out",
											  errfile=local_dir + "/mafft_log.err")
	elif context.UserFlags_dict['MSA_Software'] == 'ClustalOmega':
		run_msa_clastal(local_dir, fasta_filename, fasta_final_output)

	logger.info("FINISHED Alignment execution on %s" % fasta_filename)
	return


def run_guidance_for_seqs_and_columns(fasta_filename, guidance_work_dir1, guidance_work_dir2, fasta_final_output,
									  context, use_cache=True):
	cache_prefix = "gd_"
	dirs_for_cache = [fasta_final_output, guidance_work_dir1, guidance_work_dir2]
	logger.debug("About to check cache")

	#Check which MSA method for guidance file names MAFFT or CLUSTALO (e.g. MSA.MAFFT.Without_low_SP_Col.With_Names
	#  or   MSA.CLUSTALO.Without_low_SP_Col.With_Names)
	if context.UserFlags_dict['MSA_Software'] == 'ClustalOmega':
		guidanceMSA_name = 'CLUSTALW'
	else:
		guidanceMSA_name = 'MAFFT'

	#Michal: Guidance runs first to remove sequences and then in case seqs were removed it will run msa and guidance again for columns:
	#NEW-> Enable the OTT user to remove only Rows or only Columns:
	flag_RowCol_flag = context.UserFlags_dict['Guidance_RowCol'] # Options are: Both(default), RowsOnly, ColOnly
	if not use_cache or not use_cache_results_if_exists(fasta_filename, dirs_for_cache, context.id,
														cache_prefix):
		create_dir_if_not_exists(guidance_work_dir1)
		create_dir_if_not_exists(guidance_work_dir2)
		logger.info("Running guidance in order to remove bad sequences %s" % fasta_filename)
		run_guidance(context, fasta_filename, guidance_work_dir1)
		logger.info("DONE Running guidance in order to remove bad sequences")
		fasta_after_guidance_remove_seqs = os.path.join(guidance_work_dir1,
														"Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names")

		list_of_remove_seqs_filename = os.path.join(guidance_work_dir1, "Seqs.Orig.fas.FIXED.Removed_Seq")

		# If Guidance - Rows only: perform alignment and continue:
		if flag_RowCol_flag == 'RowsOnly':
			logger.debug("GUIDNACE filter rows only (User Flag), see filtered sequences at: %s" %list_of_remove_seqs_filename)
			logger.debug("Perform alignment on GUIDNACE output %s" %fasta_after_guidance_remove_seqs)
			last_output_from_guidance = os.path.join(guidance_work_dir1,"guidance_1_out_aligned.fasta")
			run_alignment(fasta_after_guidance_remove_seqs, last_output_from_guidance, context, guidance_work_dir1)
		elif flag_RowCol_flag == 'ColOnly':
			list_of_remove_cols_filename = os.path.join(guidance_work_dir1,'Seqs.Orig.fas.FIXED.MSA.MAFFT.Removed_Col')
			logger.debug("GUIDNACE filter columns only (USer Flag), see filtered columns at: %s" %list_of_remove_cols_filename)
			last_output_from_guidance = os.path.join(guidance_work_dir1,'MSA.MAFFT.Without_low_SP_Col.With_Names')
		else:
			#Case of both Rows and Columns:
			if os.path.exists(list_of_remove_seqs_filename) and os.path.getsize(list_of_remove_seqs_filename) == 0:
				logger.debug("No sequences were filtered using GUIDNACE - no need to run GUIDANCE again")
				last_output_from_guidance = os.path.join(guidance_work_dir1, "MSA."+guidanceMSA_name+".Without_low_SP_Col.With_Names")
			else:
				logger.info("Running guidance in order to remove bad columns %s" % fasta_after_guidance_remove_seqs)
				run_guidance(context, fasta_after_guidance_remove_seqs, guidance_work_dir2)
				fasta_after_guidance_remove_cols = os.path.join(guidance_work_dir2,
																"MSA."+guidanceMSA_name+".Without_low_SP_Col.With_Names")
				last_output_from_guidance = fasta_after_guidance_remove_cols
		# For all 3 options, copy the final filtered data file to 'fasta_final_output'
		logger.info("Copying final guidance/alignment results from %s to %s" % (last_output_from_guidance, fasta_final_output))
		shutil.copy(last_output_from_guidance, fasta_final_output)
		logger.info("DONE GUIDANCE (aligning)")
		cache_file(fasta_filename, dirs_for_cache, context.id, cache_prefix)


def get_orthomcl_table_suffix(id):
	#table_suffix = id[:14] original code - changed to enable runs with longer name to be identified (e.g. debug_runtime43 was given the name o_debug_runtime4 instead of o_debug_runtime43)
	table_suffix = id[:20]
	bad_characters = ["'", " ", "/", "\\", ":", "(", ")", "<", ">", "|", "?", "*"]
	for letter in bad_characters:
		table_suffix = table_suffix.replace(letter, "_")

	return table_suffix


def cluster_seqs(cluster_method, cluster_script_dir, id, context):
	working_dir = context.working_dir
	fasta_seq_filename = context.fasta_seq_filename
	gb_seq_filename = context.gb_seq_filename
	if (cluster_method == "orthomcl" or cluster_method == "orthomclits"):
		logger.info("Clustering sequences data using OrthoMcl ITS=" + str(context.its_support))

		if (is_cluster_results_exist(context) and not use_cached_cluster_results):
			logger.info("Cleaning old clustering results")
			clean_cluster_results(context)

		if (is_cluster_results_exist(context)):
			logger.info("Skipping clustering - clustering results found")
		else:
			logger.info("No clustering results found - starting clustering")
			create_cluster_dirs(context)

			table_suffix = get_orthomcl_table_suffix(id)
			# drop_orthomcl_tables(genus_table_suffix)
			#drop_orthomcl_tables("")
			try:
				drop_orthomcl_scehma(table_suffix)
			except Exception:
				pass

			create_orthomcl_scehma(table_suffix)
			try:
				cluster_sequences_file_with_orthomcl(fasta_seq_filename, table_suffix,
													 cluster_script_dir + "/clustering/", gb_seq_filename,
													 output_dir=context.clustering_dir, context=context)
			finally:
				drop_orthomcl_scehma(table_suffix)
	else:
		# TODO: throw exception and stop
		all_fasta_file_list = None
		final_clusters_dir = None
		logger.critical("Unknown cluster method " + cluster_method)


def is_cluster_results_exist(context):
	return os.path.exists(context.clustering_results_dir) and os.listdir(context.clustering_results_dir) != []


def clean_cluster_results(context):
	logger.info("Cleaning old cluster results")
	rm_and_create_dir(context.clustering_dir)
	rm_and_create_dir(context.cluter_its_dir)
	rm_and_create_dir(context.clustering_results_dir)


def create_cluster_dirs(context):
	logger.info("Creating clustering dirs")
	create_dir_if_not_exists(context.clustering_dir)
	create_dir_if_not_exists(context.cluter_its_dir)


def init_genus_schema(db_schema_file, db_file):
	if os.path.exists(db_file):
		logger.info("Found db file. Using existing schema %s", db_file)
	else:
		logger.info("Initializing sqlite DB file at %s using schema  %s", db_file, db_schema_file)
		conn = get_db_conn(db_file)
		curs = conn.cursor()
		with open(db_schema_file, 'r') as db_schema_handle:
			create_shema_sql = db_schema_handle.read()
			curs.executescript(create_shema_sql)

		conn.commit()
		conn.close()


def call_handle_multiple_accessions(cluster_fasta_filename, cluster_script_dir, working_dir, blast_results_filename,
									output_filename):
	if os.path.exists(output_filename):
		os.remove(output_filename)

	# Cleaning up
	rm_and_create_dir(working_dir)

	logger.info("Handling multiple accessions. rewriting %s into %s " % (cluster_fasta_filename, output_filename))
	handle_multiple_accessions(cluster_fasta_filename, blast_results_filename, working_dir, output_filename)
	# pick_one_command = "perl pickOneFromEachSpecies.pl %s %s %s %s" % (cluster_fasta_filename,blast_results_filename,working_dir,output_filename)
	#
	# logger.info("Calling the following command %s", pick_one_command)
	# exec_external_command_redirect_output(command_to_exec=pick_one_command,outfile=working_dir + "/pickOneFromEachSpecies.out",
	# errfile=working_dir + "/pickOneFromEachSpecies.err",cwd = cluster_script_dir+ "/chromCounts/miscellaneous/")
	logger.info("Finished handling multiple accessions")


def init_cluster_mult_acc_and_db(cluster_context, context, cluster_script_dir, outgrouup_selection):
	logger.info("init_cluster_mult_acc_and_db of %s" % cluster_context)
	init_cluster_mult_acc(cluster_context, context, cluster_script_dir)

	init_cluster_db(cluster_context, context, outgrouup_selection)


def init_cluster_mult_acc(cluster_context, context, cluster_script_dir):
	#
	# Handling multiple accessions
	#
	if not cluster_context.is_mult_acc_init:
		logger.info("Initializing cluster context %s " % cluster_context.cluster_id)
		logger.info("handling multiple accessions in %s" % cluster_context.all_seqs_fasta_filename)
		call_handle_multiple_accessions(cluster_context.all_seqs_fasta_filename,
										cluster_script_dir=cluster_script_dir,
										working_dir=cluster_context.multiple_accessions_work_dir,
										blast_results_filename=cluster_context.blast_results_filename,
										output_filename=cluster_context.all_seqs_fasta_filename_no_multiple_accessions)
		# handle_multiple_accessions_naiive(in_fasta=cluster_context.all_seqs_fasta_filename_no_multiple_accessions_temp, out_fasta=cluster_context.all_seqs_fasta_filename_no_multiple_accessions)
		#rm_and_create_dir(cluster_context.mrbayes_input_dir)
		cluster_context.is_mult_acc_init = True
	else:
		logger.info("cluster context %s mult acc already initialized. Skipping" % cluster_context.cluster_id)


def add_seq_to_existing_alignment(existing_msa_filename, new_seq_filename, output_filename, log_dir):
	mafft_command = "mafft --adjustdirection --addfragments %s --multipair %s > %s" % (
		existing_msa_filename, new_seq_filename, output_filename)
	retval = exec_external_command_redirect_output(mafft_command, outfile=log_dir + "/mafft-addfrag.out",
										  errfile=log_dir + "/mafft-addfrag.out")
	if retval != 0:
		raise Exception("Failed mafft cmd addfragments")

def merge_subsp_seqs_with_parent_species(context):
	# Make sure to first decide on the parent name (by removing subsp. name ext.)
	# and only then check for it's TaxID -> otherwise may cause wrong data...such as in the case of:
	# 102601|organism|Antirrhinum majus subsp. linkianum|original_taxonid|102601|original_organism|Antirrhinum linkianum
	#
	logger.info("Merging subsp. and var. sequences with their parent species sequences")
	all_seqs = list(SeqIO.parse(context.fasta_seq_filename, "fasta"))

	f_subsp_merge_data = open(context.working_dir+'/Merge_SubspVar_data.txt', 'w')
	f_subsp_merge_data.write('original_tax_id,original_name,parent_tax_id,parent_name\n')

	edited_seqs = list()
	LargeTaxon_edited_seqs = list()

	for seq in all_seqs:
		original_tax_id = getPropertyFromFastaSeqHeader(seq.description, "taxonid")
		original_name = getPropertyFromFastaSeqHeader(seq.description, "organism")

		parent_tax_id = context.parent_species_tax_id_by_subsp_tax_id[original_tax_id]
		parent_name = context.parent_species_names_by_parent_tax_id[parent_tax_id]
		f_subsp_merge_data.write(str(original_tax_id)+','+ str(original_name)+','+ str(parent_tax_id)+','+ str(parent_name)+'\n')

		restandarize_seq(seq, taxonid=parent_tax_id, original_taxonid=original_tax_id, organism=parent_name, original_organism=original_name)
		#restandarize_seq(seq, taxonid=parent_tax_id, original_subsp_taxonid=original_tax_id, organism=parent_name,original_subsp_organism=original_name)

		# For Large Taxon ids we add the sequences of organisms that are synonyms of the Accepted Large TaxID
		# to the file of the large seq.  for that we create 2 seperate lists and then append the sequences to the
		# relevant file: either the large taxId file or the general seq file
		if parent_tax_id in context.large_taxid_seqs_number_dict.keys():
			#copy the New_seq_with_AcceptedName to the LargeTaxId sequence file
			LargeTaxon_edited_seqs.append(seq)
		else:
			#copy to the regular sequence file
			edited_seqs.append(seq)
		#edited_seqs.append(seq)

	shutil.copyfile(context.fasta_seq_filename, context.fasta_seq_merged_subsp_filename)
	write_standarize_fasta(context.fasta_seq_filename, edited_seqs)
	if LargeTaxon_edited_seqs:	# if there's taxa to add to the large taxon (syn of Large taxa:
		write_standarize_fasta(context.largeSynToAdd_file, LargeTaxon_edited_seqs)
		append_seq_to_largeTaxId(context)
	logger.info("Finished merging subsp. and var. sequences with their parent species sequences")

#-------------------------------------------------------------------------------------------------------------------
#This function will merge the taxons which belong to the same Accepted name: so each synonym of an accepted name
#will be replaced by the taxID of the accepted and also the name. They will be considered as one
def merge_syn_seqs_with_accepted_species(context):
	logger.info("Merging synonyms sequences with their accepted species sequences")
	all_seqs = list(SeqIO.parse(context.fasta_seq_filename, "fasta"))
	# Using the following dict for syn->acc transfer
	# context.syn_acc_dict - taxid vs name
	# context.syn_acc_TaxID_dict - taxid vs taxid
	##(syn_acc_dict): e.g.          '309346': 'Linaria genistifolia',
	##(syn_acc_TaxID_dict): e.g.    '1851693': '309346',


	#Open file to save all data for subsp merge:
	f_SynAcc_merge_data = open(context.working_dir+'/Merge_SynonymToAccepted_data.txt', 'w')
	f_SynAcc_merge_data.write('original_tax_id,original_name,parent_tax_id,parent_name\n')

	edited_seqs = list()
	LargeTaxon_edited_seqs = list()

	already_in_list=[]
	#Create a file that will summarize the input to the pipeline:
	#Spc_data_f = open(context.working_dir+'/Summary_Spc_with_MergeData.txt','w')
	for seq in all_seqs:
		original_tax_id = getPropertyFromFastaSeqHeader(seq.description, "taxonid")
		original_name = getPropertyFromFastaSeqHeader(seq.description, "organism")

		#Need to check if it's a synonym and if the Accepted name has a TaxID:
		#First check if Synonym:
		#Accepted_name,AcceptedId = get_Accepted_name(original_tax_id)
		logger.info("Check for replaceing synonym, with accepted: original_tax_id %s original_name %s" % (original_tax_id,original_name))
		if original_tax_id in context.syn_acc_TaxID_dict.keys():
			parent_tax_id = context.syn_acc_TaxID_dict[original_tax_id]
			parent_name= context.syn_acc_dict[parent_tax_id]
			logger.info("%s(%s) is a synonym of %s(%s). Replacing sequence description." % (original_name,original_tax_id, parent_name,parent_tax_id))
			logger.debug("Original Name %s, Accepted Name %s" %(original_name,parent_name))
			if str(original_name) not in already_in_list:
				already_in_list.append(str(original_name))
				f_SynAcc_merge_data.write(str(original_tax_id)+','+ str(original_name)+','+ str(parent_tax_id)+','+ str(parent_name)+'\n')
				context.original_vs_newName_dict[str(original_name)]=str(parent_name)
				context.spc_tax_afterMatch_dict[str(parent_name)]=str(parent_tax_id)
		else:
			parent_tax_id = original_tax_id
			parent_name = original_name
			context.original_vs_newName_dict[str(original_name)]=str(parent_name)
			context.spc_tax_afterMatch_dict[str(parent_name)]=str(parent_tax_id)

		restandarize_seq(seq, taxonid=parent_tax_id, original_taxonid=original_tax_id, organism=parent_name, original_organism=original_name)

		# For Large Taxon ids we add the sequences of organisms that are synonyms of the Accepted Large TaxID
		# to the file of the large seq.  for that we create 2 seperate lists and then append the sequences to the
		# relevant file: either the large taxId file or the general seq file
		if parent_tax_id in context.large_taxid_seqs_number_dict.keys():
			#copy the New_seq_with_AcceptedName to the LargeTaxId sequence file
			LargeTaxon_edited_seqs.append(seq)
		else:
			#copy to the regular sequence file
			edited_seqs.append(seq)


	shutil.copyfile(context.fasta_seq_filename, context.fasta_seq_merged_syn_filename)
	write_standarize_fasta(context.fasta_seq_filename, edited_seqs)
	if LargeTaxon_edited_seqs:	# if there's taxa to add to the large taxon (syn of Large taxa:
		write_standarize_fasta(context.largeSynToAdd_file, LargeTaxon_edited_seqs)
		append_seq_to_largeTaxId(context)

	#Now we need to handle all Large species which are synonyms of regular species:
	large_renewd_list = []
	for large_taxa in context.large_taxid_seqs_number_dict.keys():
		large_fasta_file = context.largeSpeciesDir + '/seqs_for_LargeTaxId_' + str(large_taxa) + '.fasta'
		curr_large_seqs = list(SeqIO.parse(large_fasta_file, "fasta"))
		for seq in curr_large_seqs:
			if large_taxa in context.syn_acc_TaxID_dict.keys():
				parent_tax_id = context.syn_acc_TaxID_dict[large_taxa]
				parent_name= context.syn_acc_dict[parent_tax_id]
				logger.info("(Large) %s(%s) is a synonym of %s(%s). Replacing sequence description." % (original_name,original_tax_id, parent_name,parent_tax_id))
				logger.debug("(Large) Original Name %s, Accepted Name %s" %(original_name,parent_name))
				if str(original_name) not in already_in_list:
					already_in_list.append(str(original_name))
					f_SynAcc_merge_data.write(str(original_tax_id)+','+ str(original_name)+','+ str(parent_tax_id)+','+ str(parent_name)+'\n')
					context.original_vs_newName_dict[str(original_name)]=str(parent_name)
					context.spc_tax_afterMatch_dict[str(parent_name)]=str(parent_tax_id)
			else:
				parent_tax_id = original_tax_id
				parent_name = original_name
				context.original_vs_newName_dict[str(original_name)]=str(parent_name)
				context.spc_tax_afterMatch_dict[str(parent_name)]=str(parent_tax_id)
			#rewrite the seq with the correct accepted name:
			restandarize_seq(seq, taxonid=parent_tax_id, original_taxonid=original_tax_id, organism=parent_name, original_organism=original_name)
			large_renewd_list.append(seq)
		if large_renewd_list:
			shutil.copyfile(large_fasta_file, large_fasta_file+'_origin')
			write_standarize_fasta(large_fasta_file, large_renewd_list)


	logger.info("Finished merging synonyms sequences with their accepted species sequences")
	#sys.exit()
	return

def perform_filter_msa(ploidb_context, msa_file):
	#Run Filter method to avoid large msa files:
	filter_method = ploidb_context.UserFlags_dict['FilterMSA_Method']
	if filter_method == 'Trimal':
		Trimal_cf = ploidb_context.UserFlags_dict['Trimal_CutOff']
		runTrimal(msa_file, Trimal_cf)
		logger.info("Files after Trimal at: %s" % msa_file)
	elif filter_method == 'GUIDANCE':
		run_guidance_for_seqs_and_columns(msa_file, cluster_context.guidance_concat_work_dir1,
										  cluster_context.guidance_concat_work_dir2,
										  cluster_context.all_seqs_concat_aligned_fasta_filename, context)
	elif filter_method == 'gBlocks':
		runGblocks(msa_file, ploidb_context.UserFlags_dict)
	else:
		logger.debug("User didn't select any MSA filter method")
		return
	return

def perform_filter_ITS_msa(ploidb_context, msa_file):
	#Run Filter method to avoid large msa files:
	filter_method = ploidb_context.UserFlags_dict['FilterMSA_Method']
	if filter_method == 'Trimal':
		Trimal_cf = ploidb_context.UserFlags_dict['Trimal_CutOff']
		runTrimal(msa_file, Trimal_cf)
		logger.info("Files after Trimal at: %s" % msa_file)
	elif filter_method == 'GUIDANCE':
		guidance_concat_work_dir1 = ploidb_context.cluter_its_dir + '/guidance1'
		guidance_concat_work_dir2 = ploidb_context.cluter_its_dir + '/guidance2'
		create_dir_if_not_exists(guidance_concat_work_dir1)
		create_dir_if_not_exists(guidance_concat_work_dir2)
		run_guidance_for_seqs_and_columns(msa_file, guidance_concat_work_dir1,
										  guidance_concat_work_dir2,
										  cluster_context.all_seqs_concat_aligned_fasta_filename, context)
	elif filter_method == 'gBlocks':
		runGblocks(msa_file, ploidb_context.UserFlags_dict)
	else:
		logger.debug("User didn't select any MSA filter method")
		return
	return


#New outgroup flow: Michal 2017
def process_seq_and_create_concat_file(context, genus_name, outgrouup_selection):
	if outgrouup_selection == 'None':
		logger.info("No Outgroup selection for  " + genus_name)
		init_NO_common_outgroup(context, outgrouup_selection)
	#elif context.UserFlags_dict['Outgroup_User'] != 'None':
	elif 'Outgroup_Flag' in context.UserFlags_dict:
		if context.UserFlags_dict['Outgroup_Flag'] == 'User':
			logger.info("Outgroup selected by the user: %s " % context.UserFlags_dict['Outgroup_User'])
			init_NO_common_outgroup(context, context.UserOutgroupName)
		elif context.UserFlags_dict['Outgroup_Flag'] == 'Single':
			logger.info("Initializing common outgroup  for  " + genus_name)
			init_common_outgroup(context)

	# For Outgroup option, check if seqs are ready:
	if is_seqs_for_concat_tree_exists(context):
		logger.info("All sequences for concat tree creating exist. Skipping their creation")
	else:
		logger.info("Writing seq files that contains the concat outgroup")
		for idx, cluster_context in enumerate(context.cluster_contexts_for_concat_tree):
			#ITS cluster only:
			if cluster_context.is_its_cluster:
				if outgrouup_selection == 'None':
					shutil.copyfile(cluster_context.all_seqs_fasta_filename_no_multiple_accessions,
									cluster_context.all_seqs_concat_aligned_fasta_filename)
				else:
					#in case of an outgroup, add it using mafft:
					if cluster_context.no_of_outgroup_sequence_for_concat > 0:
						logger.info("Calling mafft in order to add ITS outgroup to the concat sequence fasta")
						add_seq_to_existing_alignment(cluster_context.all_seqs_fasta_filename_no_multiple_accessions,
													  cluster_context.outgroup_sequence_for_concat_filename,
													  cluster_context.all_seqs_with_concat_outgroup_filename,
													  cluster_context.cluster_concat_work_dir)
						shutil.copyfile(cluster_context.all_seqs_with_concat_outgroup_filename,
										cluster_context.all_seqs_concat_aligned_fasta_filename)
					else:
						logger.info("No outgroup - copying %s to %s" % (
							cluster_context.all_seqs_fasta_filename_no_multiple_accessions,
							cluster_context.all_seqs_concat_aligned_fasta_filename))
						shutil.copyfile(cluster_context.all_seqs_fasta_filename_no_multiple_accessions,
										cluster_context.all_seqs_concat_aligned_fasta_filename)
					logger.info("DONE Adding ITS outgroup")
			else:
				#Handling other clusters, part of ITS:
				logger.info(
					"Writing seq file that contains the concat outgroup for into %s" % cluster_context.all_seqs_with_concat_outgroup_filename)
				if outgrouup_selection == 'None':
					shutil.copyfile(cluster_context.all_seqs_fasta_filename_no_multiple_accessions,
									cluster_context.all_seqs_with_concat_outgroup_filename)
				else:
					# Clusters with outgroup:
					if cluster_context.no_of_outgroup_sequence_for_concat > 0:
						merge_seq_files([cluster_context.all_seqs_fasta_filename_no_multiple_accessions,
										 cluster_context.outgroup_sequence_for_concat_filename],
										cluster_context.all_seqs_with_concat_outgroup_filename, 'False')
					else:
						shutil.copyfile(cluster_context.all_seqs_fasta_filename_no_multiple_accessions,
										cluster_context.all_seqs_with_concat_outgroup_filename)
				logger.info("START handling reversed seqs for concat")
				adjust_direction_fasta_using_mafft([cluster_context.all_seqs_with_concat_outgroup_filename])
				logger.info("DONE handling reversed seqs for concat")
				fasta_for_alignment = cluster_context.all_seqs_with_concat_outgroup_filename

				if os.path.exists(cluster_context.all_seqs_concat_aligned_fasta_filename):
					logger.debug(
						"File %s already exists - skipping Alignment" % cluster_context.all_seqs_concat_aligned_fasta_filename)
				else:
					# Perform Alignment:
					if context.UserFlags_dict['FilterMSA_Method'] == 'GUIDANCE':
						run_guidance_for_seqs_and_columns(fasta_for_alignment,
														  cluster_context.guidance_concat_work_dir1,
														  cluster_context.guidance_concat_work_dir2,
														  cluster_context.all_seqs_concat_aligned_fasta_filename,
														  context)
					else:
						run_alignment(fasta_for_alignment, cluster_context.all_seqs_concat_aligned_fasta_filename,
									  context, cluster_context.cluster_concat_work_dir)
						if context.UserFlags_dict['FilterMSA_Method'] is not 'None':
							# Check MSA Filter Methods:
							logger.debug(
								"User selected Filter MSA Method is: %s" % context.UserFlags_dict['FilterMSA_Method'])
							perform_filter_msa(context, cluster_context.all_seqs_concat_aligned_fasta_filename)

			logger.info("Rewriting %s to %s - using organism name as seq description" % (
				cluster_context.all_seqs_concat_aligned_fasta_filename,
				cluster_context.all_seqs_for_concat_tree_creation_fasta_filename))
			names_dict1 = rewrite_fasta_with_organism_name_as_desc(
				cluster_context.all_seqs_concat_aligned_fasta_filename,
				cluster_context.all_seqs_for_concat_tree_creation_fasta_filename, False)
			names_dict2 = dict()
			if cluster_context.no_of_outgroup_sequence_for_concat > 0:
				names_dict2 = rewrite_fasta_with_organism_name_as_desc(
					cluster_context.outgroup_sequence_for_concat_filename,
					cluster_context.outgroup_sequence_for_concat_short_desc_filename, False)
			else:
				logger.debug("No outgroup found for cluster %s" % cluster_context)

			all_names_dict = names_dict1.copy()
			all_names_dict.update(names_dict2)
			with open(context.rewrite_names_dict_filename, 'w') as names_handle:
				json.dump(all_names_dict, names_handle)

			logger.info("DONE Rewriting")

		logger.info("DONE Writing seq files that contains the concat outgroup")

	context.init_chloroplast_cluster()
	context.init_mitochondrial_cluster()
	context.init_nuc_cluster()

	logger.info("Concatenating alignments for  " + genus_name)
	concat_alignment_SingleOutgroup(context)

	logger.info("DONE Concatenating alignments for  " + genus_name)
	#Create_species_list included in the Alignment:
	with open(context.concat_seqs_fasta_filename, 'r') as f_alignment:
		f_finalNames = open(context.final_species_file, 'w')
		if outgrouup_selection == 'None':
			for line in f_alignment:
				if '>' in line:
					line = line.strip()
					spc_name = line.replace('>', '').replace('_', ' ')
					f_finalNames.write(spc_name + '\n')
					context.list_final_species_names.append(spc_name)
		else:
			logger.debug("Creating Final list from concat file:")
			for line in f_alignment:
				if '>' in line and context.outgroupSelection not in line:
					line = line.strip()
					logger.debug(line)
					spc_name = line.replace('>', '').replace('_', ' ')
					#f_finalNames.write('%s,%s\n' %(context.Pipe_input_Spc_Names_vs_TaxId[spc_name],spc_name))
					#f_finalNames.write('%s,%s\n' %(context.species_names_by_tax_id[spc_name],spc_name))
					f_finalNames.write('%s\n' %(spc_name))
					context.list_final_species_names.append(spc_name)
		f_finalNames.close()
		f_alignment.close()


# This is the main method - it invokes the rest of the methods in the pipeline
def verify_cluster_data_is_enough(context, init_check=False):
	is_ok = True
	reason = None

	path_helper = PathHelper(context.working_dir, context.id)
	# if os.path.exists(path_helper.orthomcl_mcl_input):
	#	mcl_input_size = os.path.getsize(path_helper.orthomcl_mcl_input)
	#	if mcl_input_size == 0:
	#		is_ok = False
	#		reason = "OrthoMCL input is empty - not enough data after clustering"

	min_species_in_cluster = int(ott_config['general']['min_species_in_cluster'])
	unresolved_species_no = 0
	if context.unresolved_species_names is not None:
		unresolved_species_no = len(context.unresolved_species_names)

	if init_check:
		different_species_after_name_resolved = get_number_of_different_species(context.fasta_seq_filename)
		if different_species_after_name_resolved < min_species_in_cluster:
			is_ok = False
			reason = "No data after name resolution."

	else:
		if len(context.cluster_contexts) == 0:
			is_ok = False
			reason = "No clusters were found to work with."
		elif context.cluster_contexts[0].get_number_of_different_species() < min_species_in_cluster:
			is_ok = False
			reason = "No clusters with more than %d species were found." % min_species_in_cluster

	if reason is not None:
		reason += " Removed by name resolution %d" % unresolved_species_no
	if not is_ok and reason is not None:
		with open(path_helper.pre_tree_stopped, "w") as handle:
			handle.write(reason)

	return is_ok


def create_summary_zip(context):
	#Copy user input file to summary dir:
	if os.path.exists(context.working_dir+'/userInput.txt'):
		shutil.copyfile(context.working_dir+'/userInput.txt',context.summary_dir+'/userInput.txt')
	#Create zipped file of summary dir
	summary_zipName = context.working_dir + '/OneTwoTree_Output_' + context.id + '.zip'
	logger.debug("Zip Summary dir: %s to: %s" % (context.summary_dir, summary_zipName))
	sum_zip = zipfile.ZipFile(summary_zipName, mode='w')
	for root, dirs, files in os.walk(context.summary_dir):
		for file in files:
			sum_zip.write(os.path.join(root, file), file)
	sum_zip.close()
	return


def create_LocusTree_diffMethods(cluster_id, Tree_Method, context):
	alg_Locus_file = context.working_dir + '/concat/' + cluster_id + '/seqs-organism-concat.fasta'
	if Tree_Method == "ExaML":
		logger.debug("# ----------> Tree Method chosen: ExaML")
		xml_model = context.UserFlags_dict['User_ExaML_Model']
		Build_examl_LocusTree(alg_Locus_file, context, xml_model, cluster_id)
	elif Tree_Method == "RAxML":
		logger.debug("# ----------> Tree Method chosen: RAxML")
		xml_model = context.UserFlags_dict['User_RAxML_Model']
		Build_raxml_LocusTree(alg_Locus_file, context, xml_model, cluster_id)
	else:
		logger.debug("# ----------> ERROR Tree Method doesn't exist: %s" % Tree_Method)
		statusFail_LogFile(context, "Tree Method doesn't exist")
		return
	create_summary_zip(context)
	#Completed the flow - PAssed all
	with open(context.final_status, 'w') as f_status:
		f_status.write("PASS - Tree Created Successfully")
		f_pass = open(context.working_dir+'/JOB_PASS.txt','w')
		f_pass.close()
	return


def create_Tree_diffMethods(f_stat_check, Tree_Method, context):
	#Build_MB_Tree(f_end_status_name,context)
	outgorup_flag = context.UserFlags_dict['Outgroup_Flag']


	NameDateFlag='False'
	if "UltraMetricFlag" in context.UserFlags_dict:
		UltrametricFlag = context.UserFlags_dict["UltraMetricFlag"]
	else:
		UltrametricFlag = 'No'
	TreeUserchoise = 0
	mb_nodeDate_txt = 'Empty'
	NameDate_missingName_Flag = 'Yes'
	#Check for Node dating params:
	if os.path.exists(context.NodeDatingFile):
		f_NodeDate = open(context.NodeDatingFile, 'r')
		if is_empty_file(f_NodeDate) != 0:
		#if os.stat(context.NodeDatingFile).st_size != 0:
			logger.debug("User gave Node dating splits data...reading file:")
			first_name_list=[]
			second_name_list=[]
			min_date_list=[]
			max_date_list=[]
			NumOfNodeDates=0
			#with open(context.NodeDatingFile,'r') as f_NodeDate:
			for line in f_NodeDate:
				if line.strip():
					new_line=list([x.strip() for x in line.split(',')])
					logger.debug(new_line)
					first_name_list.append(new_line[0])
					second_name_list.append(new_line[1])
					min_date_list.append(new_line[2])
					max_date_list.append(new_line[3])
					NumOfNodeDates+=1
			logger.debug("Number of splits entered by the user are: %d" %NumOfNodeDates)
			msa_conct_length = return_msa_length(context.concat_seqs_report_filename)
			#check if names are on the tree:
			NameDate_missingName_Flag='No'
			idx=0
			while idx < NumOfNodeDates:
				if first_name_list[idx] in context.list_final_species_names and second_name_list[idx] in context.list_final_species_names:
					logger.debug("%d Names for node dating: %s,%s were found on the final tree" %(idx,first_name_list[idx],second_name_list[idx]))
					idx+=1
				else:
					logger.debug("%d Names for node dating: %s,%s were NOT found on the final tree, Node dating tree will not be generated" %(idx,first_name_list[idx],second_name_list[idx]))
					NameDate_missingName_Flag='Yes'
					idx += 1
			NameDateFlag='True'
			create_FileNameDate_config(context,msa_conct_length,NumOfNodeDates,first_name_list,second_name_list,min_date_list,max_date_list)
			idx=0
			while idx < NumOfNodeDates:
				logger.debug("%s,%s,%s,%s" %(first_name_list[idx],second_name_list[idx],min_date_list[idx],max_date_list[idx]))
				idx+=1

			#Check if names are on the final msa and if not no need to create the data for the config
			if NameDate_missingName_Flag == 'No':
				mb_nodeDate_txt = create_NodeDate_paragraph_MB(NumOfNodeDates,first_name_list,second_name_list,min_date_list,max_date_list)
			if NumOfNodeDates == 0:
				mb_nodeDate_txt = 'Empty'


	##------------------------------------------- MrBayes -----------------------------------------------
	if Tree_Method == "MrBayes":
		if context.UserFlags_dict["Constraint_Flag"] == 'On':   # User would like to add a constraint tree so we need to create the blk to add to the mb config file
			create_mb_blk_constraint(context)
		logger.debug("# ----------> Tree Method chosen: MrBayes")
		if Build_MB_Tree(f_stat_check, context, mb_nodeDate_txt) == 'fail':
			return
	##------------------------------------------- EXAML -----------------------------------------------
	elif Tree_Method == "ExaML":
		logger.debug("# ----------> Tree Method chosen: ExaML")
		xml_model = context.UserFlags_dict['User_ExaML_Model']
		if context.UserFlags_dict['Use_BootS_ExOn'] == 'Off':
			Build_examl_Tree(context.concat_seqs_fasta_filename, context, xml_model)
		else:
			bootstrap_val = str(context.UserFlags_dict['Examl_BootS_vals'])
			Build_examl_Boots_Tree(context.concat_seqs_fasta_filename, context, xml_model,bootstrap_val,outgorup_flag)

		#Ultrametric without Node Dating option
		#retval=0
		#Ultra / Constriant etc
		finalTree = context.summary_dir + '/Result_Tree_' + context.id + '.tre'
		if "NodeDateMethod" in context.UserFlags_dict:
			ml_tree_method = context.UserFlags_dict["NodeDateMethod"]
		if UltrametricFlag == 'Yes' and NameDateFlag == 'False' and os.path.exists(finalTree) and NameDate_missingName_Flag=='Yes':
			perform_ML_UltrametricNoDate(context,ml_tree_method,finalTree)
		#Node Dating option
		elif UltrametricFlag == 'Yes' and NameDateFlag == 'True' and os.path.exists(finalTree) and NameDate_missingName_Flag=='No':
			perform_ML_UltrametricWithDate(context,ml_tree_method,finalTree,NameDate_missingName_Flag)
			#perform_ML_UltrametricWithDate(context,ml_tree_method,tree_file_path,NameDate_missingName_Flag)


	##------------------------------------------- RAXML -----------------------------------------------
	elif Tree_Method == "RAxML":
		logger.debug("# ----------> Tree Method chosen: RAxML")
		xml_model = context.UserFlags_dict['User_RAxML_Model']
		if context.UserFlags_dict['Use_BootS_RxOn'] == 'Off':
			Build_raxml_Tree(context.concat_seqs_fasta_filename, context, xml_model)
		else:
			bootstrap_val = str(context.UserFlags_dict['Raxml_BootS_vals'])
			Build_raxml_Boots_Tree(context.concat_seqs_fasta_filename, context, xml_model,bootstrap_val,outgorup_flag)
		#Ultra / Constriant etc
		if "NodeDateMethod" in context.UserFlags_dict:
			ml_tree_method = context.UserFlags_dict["NodeDateMethod"]
		finalTree = context.summary_dir + '/Result_Tree_' + context.id + '.tre'
		logger.debug("UltrametricFlag - %s, NameDateFlag - %s,NameDate_missingName_Flag - %s, finalTree - %s" %(UltrametricFlag,NameDateFlag,NameDate_missingName_Flag,finalTree))
		if UltrametricFlag == 'Yes' and NameDateFlag == 'False' and os.path.exists(finalTree) and NameDate_missingName_Flag=='Yes':
			perform_ML_UltrametricNoDate(context,ml_tree_method,finalTree)
		#Node Dating option
		elif UltrametricFlag == 'Yes' and NameDateFlag == 'True' and os.path.exists(finalTree) and NameDate_missingName_Flag=='No':
			perform_ML_UltrametricWithDate(context,ml_tree_method,finalTree,NameDate_missingName_Flag)


	##------------------------------------------- ERROR -----------------------------------------------
	else:
		logger.debug("# ----------> ERROR Tree Method doesn't exist: %s" % Tree_Method)
		statusFail_LogFile(context, "Tree Method doesn't exist")
		return
	#Zip_summary dir:
	if os.path.exists(context.summary_dir + '/Result_Tree_' + context.id + '.tre'):
		final_tree = context.summary_dir + '/Result_Tree_' + context.id + '.tre'
		final_tree_pdf = context.summary_dir + '/Result_Tree_' + context.id + '.pdf'
		scripts_dir =  ott_config['general']['OTT_MAIN'] + 'ott_scripts/'
		r_path = ott_config['diff_soft']['R_PATH']
		createPdfTree_Rcmd = (r_path + " CMD BATCH '--args working_dir=" + '"' + context.summary_dir + '"' + " tree=" + '"'
		+ final_tree + '"' + " output=" + '"' + final_tree_pdf + '"' + "' " + scripts_dir + "plot_pdf.R")
		logger.debug("R script to create pdf file: %s" %createPdfTree_Rcmd)
		exec_external_command_redirect_output(createPdfTree_Rcmd)
		create_summary_zip(context)
		#Completed the flow - Passed all
		with open(context.final_status, 'w') as f_status:
			f_status.write("PASS - Tree Created Successfully")
			f_pass = open(context.working_dir+'/JOB_PASS.txt','w')
			f_pass.close()
	else:
		with open(context.final_status, 'w') as f_status:
			f_status.write("! Tree was not created, please send us your job number for more details")

	return


def statusFail_LogFile(context, status_line):
	logger.info(
		"=======================================================================================================")
	logger.info("===                 %s	         " % status_line)
	logger.info(
		"=======================================================================================================")
	#Result indication file:
	with open(context.final_status, 'a') as f_status:
		f_status.write("Failed - %s \n" % status_line)
	return


# ------------------------------------------------------------------------------------------------------------------------------------------------------------
#
#                                                               Main function of Pipeline:
#
# ------------------------------------------------------------------------------------------------------------------------------------------------------------

def generate_taxa_tree(id, taxa_list, working_dir, log_filename, debug_filename, cluster_method=None,
					   should_ignore_subsp=False, should_merge_subsp=False):
	# init log
	if log_filename is None:
		log_filename = working_dir + "/" + id + ".log"
	if debug_filename is None:
		debug_filename = working_dir + "/" + id + "-debug.log"
	error_filename = working_dir + "/" + id + "-error.log"

	initLogger(log_filename, debug_filename, error_filename)
	initCommandLineLogger()


	logger.info(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Start processing job %s >>>>>>>>>>>>>>>>>>>>>", id)
	p = PathHelper(working_dir, id)
	open(p.pre_tree_started, 'a').close()

	if os.path.exists(p.pre_tree_stopped):
		os.remove(p.pre_tree_stopped)

	init_ott_config(args.config_filename)
	logger.info("ott_config=%s" % ott_config)
	cluster_script_dir = ott_config['general']['OTT_MAIN']
	#Check if to diable debug printouts:
	if ott_config['general']['DEBUG_PRINT'] == 'off':
		logging.disable(logging.DEBUG)
	context = PloiDbContext(id=id, working_dir=working_dir, cluster_script_dir=cluster_script_dir,
							cluster_method=cluster_method, taxa_list=taxa_list)

	# Update User Flags from file - params
	create_dir_if_not_exists(context.summary_dir)
	create_dir_if_not_exists(context.debug_dir)
	UpdateContextFlags(context)
	context.should_run_guidance = context.UserFlags_dict['FilterMSA_Method']
	open(context.summary_file, 'w')
	outgrouup_selection = context.UserFlags_dict['Outgroup_Flag']
	context.UserFlag_MSAmethod = context.UserFlags_dict['MSA_Software']
	context.UserFlag_mbUserModel = context.UserFlags_dict['User_MrBayes_Model']
	context.ott_scripts_dir = ott_config['general']['OTT_MAIN'] + '/ott_scripts/'
	#Stand alone version :
	context.standAlone_Flag = ott_config['general']['SA_VERSION']
	#Data validation:
	if 'Outgroup_User' in context.UserFlags_dict:
		if context.UserFlags_dict['Outgroup_Flag'] == 'User':
			if not context.UserFlags_dict['Outgroup_User']:	# User forgot to insert outgroup name
				status_line = 'missing outgroup name, resend and make sure to fill in the user defined outgroup name'
				logger.debug(status_line)
				with open(context.final_status, 'a') as f_status:
					f_status.write("Failed - %s \n" % status_line)
				raise Exception("Failed - missing user defined value")


	# Diff running params:
	context.cluster_script_dir = cluster_script_dir
	context.id = id
	context.working_dir = os.path.abspath(working_dir)
	context.should_ignore_subsp = should_ignore_subsp
	context.should_merge_subsp = should_merge_subsp
	db_schema_filename = context.ott_scripts_dir + "/genusSchema.sql"

	# Creating ploidb.db file:
	init_genus_schema(db_schema_filename, context.db_file)
	init_genus_schema(db_schema_filename, context.ploidb_db_file)
	context.start_time = datetime.datetime.now()

	#Validate input fies:
	if context.UserFlags_dict["Constraint_Flag"] == "On":
		constraintTree_file = context.working_dir+'/ConstraintTree_user.txt'
		try:
			tree = Phylo.read(constraintTree_file, 'newick')
		except:
			logger.debug("Constraint Tree format must be Newick, please check it and try again")
			statusFail_LogFile(context, 'Constraint Tree format must be Newick, please check it and try again')
			return


	#check file to be reg txt file:
	for item in context.taxa_list:
		if item.startswith(">"):
			statusFail_LogFile(context, 'Your input file is in the wrong format, should be a txt file with a list of taxa names/TaxIds')
			with open(context.working_dir+'/calc_time_vars.txt','w') as f_time_calc:
				f_time_calc.write(' ')
				f_time_calc.close()
			return


	# Check if Rerun
	#-------------------------------------

	context.RerunID = context.UserFlags_dict['OriginJobID']

	if context.RerunID == 'None':
		logger.info("Rerun Job ID is %s....continue regular OTT" % context.RerunID)
	else:
		logger.info("Rerun on Job ID: %s" % context.RerunID)
		Rerun_flow(context.RerunID, context)
		logger.info("Done with Rerun flow !!!!!!!!!!!!!!!!!")
		return

	# Check if Clustering option - OTT site option: runningOptions:clustering
	#-------------------------------------
	context.OTTOptions = context.UserFlags_dict['runningOptions']
	logger.debug("User running option is: %s" % context.OTTOptions)
	# 'phylogeny reconstruction' or 'clustering'


	# ------------------------------------
	#Check if Alignment is done - continue to TreeGen:
	f_stat_check = context.debug_dir + '/AlignEnd_FinishedExecutionWithConcatFile'

	#try:
	if os.path.exists(f_stat_check):
		context.outgroupSelection = rerun_get_selected_outgroup(context.working_dir + "/OutgroupSelection.txt")
		#Copy alignment file to Summary dir:
		shutil.copyfile(context.concat_seqs_fasta_filename,
						context.summary_dir + '/' + context.id + "-concat-aligned.fasta")
		logger.info(
			"======================     Alignment DONE - Start Tree Generation	====================================")
		Create_Tree(context,f_stat_check)
		return
	else:
		logger.info(
			"======================     START PERFORMING ALIGNMENT SECTION	====================================")

	# In case of User outgroup, add it to taxa list:
	db_name = ott_config['general']['DB_DIR'] + ott_config['concat_headir']['NCBI_NAMES_NR']
	if 'Outgroup_User' in context.UserFlags_dict:
		if context.UserFlags_dict['Outgroup_Flag'] == 'User':
			context.taxa_list.append(context.UserFlags_dict['Outgroup_User'])
			logger.debug("Added User outgroup: %s to taxa list" % context.UserFlags_dict['Outgroup_User'])
			User_outgroup = context.UserFlags_dict['Outgroup_User']
			if User_outgroup.isdigit():
				taxon_item = User_outgroup.rstrip()
				query_ncbi = "SELECT Name,TaxID from ncbi_names_db where TaxID like '%s' and Type like 'scientific name'" % taxon_item
			else:
				taxon_item = User_outgroup.rstrip()
				query_ncbi = "SELECT Name,TaxID from ncbi_names_db where Name like '%s' and Type like 'scientific name'" % taxon_item
			ncbi_res = query_ncbi_db(query_ncbi, db_name)
			if ncbi_res:
				for line in ncbi_res:
					Name = str(line[0])
					TaxId = str(line[1])
					break
				context.UserOutgroupTaxId = TaxId
				context.UserOutgroupName = Name
				logger.debug("User outgroup: %s, TaxId %s" %(Name,TaxId))
				logger.debug(Name)
				context.outgroupSelection = Name.replace(' ','_')
				context.UserOutgroupTaxId = TaxId
			else:
				logger.debug("User Outgroup was not found in NCBI genabnk database")
				status_line = 'User Outgroup was not found in NCBI genabnk database'
				logger.debug(status_line)
				with open(context.final_status, 'a') as f_status:
					f_status.write("Failed - %s \n" % status_line)
				raise Exception("Failed - missing user defined value")


	#New species list - NR - Still under DEBUG:
	#---------------------------------------------------------------------------------------------------
	# 1. None - check if names are in DB -> Yes, add to names list, No, Add to NOT FOUND file
	# 2. Spelling - check name spelling according to chosen DB use OTT_NR.py
	# 3. Matched - prepare dict Accepted for each synonym

	# Step 1:
	return_species_list_debug(context.taxa_list,context)

	#Construct species list and taxIds:
	logger.info("species List:")
	logger.info(context.species_list_names)
	logger.info("species Ids:")
	logger.info(context.species_list_ids)
	if context.UserFlags_dict['NameResType'] == 'SpellingCorrection' or context.UserFlags_dict['NameResType'] == 'MatchedNames':
		logger.info("For each TaxId, Accepted name (syn_acc_dict)")
		logger.info(context.syn_acc_dict)
		logger.info("For each TaxId, Accepted name (syn_acc_TaxID_dict)")
		logger.info(context.syn_acc_TaxID_dict)

	if len(context.species_list_ids) < 1:
		statusFail_LogFile(context, 'None of your input taxa was found in Genbank, Please check your input list')
		with open(context.working_dir+'/calc_time_vars.txt','w') as f_time_calc:
			f_time_calc.write(' ')
			f_time_calc.close()
		return
	if len(context.species_list_ids) < 5:
		statusFail_LogFile(context, 'Less than 5 species found')
		with open(context.working_dir+'/calc_time_vars.txt','w') as f_time_calc:
			f_time_calc.write(' ')
			f_time_calc.close()
		return

	# ------------------------------------
	# Checking if species list length is not to long:
	max_number_of_species = int(ott_config['general']['max_number_of_species'])
	species_list_length = len(context.species_list_ids)
	logger.debug(
		"Checking species list length. Species list contains %d species, max number of species allowed is %d" % (
			species_list_length, max_number_of_species))
	if species_list_length > max_number_of_species:
		with open(context.summary_file, 'a') as f_sum:
			f_sum.write("EndStatus: Species list is Too long: %d species found, max number allowed is %d. Please contact us to enable this run.\n" %
						(species_list_length, max_number_of_species))
		raise Exception("Species list is to long. It contains %d species, max number of species allowed is %d" % (
			species_list_length, max_number_of_species))

	# ------------------------------------
	# Getting sequences for all taxa list taxon_names
	shutil.copyfile(context.species_names_taxid_file, context.summary_dir + '/species_list_ids.txt')
	# The seqs are saved in .gb and .fasta files under the working dir
	if os.path.exists(context.fasta_seq_filename):
		logger.info("## Using existing sequences fasta file: %s !!!" % context.fasta_seq_filename)
	else:
		get_taxa_sequences(context, context.species_list_ids, fasta_filename=context.fasta_seq_filename,
						   gb_filename=context.gb_seq_filename, p=p)
	logger.info("DONE Getting sequences data for taxa")

	#Update file with all TaxIds without Data in Genbank (wither they don't exist or they are of higher rank type:
	with open(context.TaxId_NoData_list_file,'a') as f_NoData:
		for TaxId in context.taxIds_not_in_genbank_list:
			f_NoData.write(TaxId + '\n')

	#MD - update context.species_list_names with the species participating in the alignment
	all_seqs = list(SeqIO.parse(context.fasta_seq_filename, "fasta"))
	species_names_from_fasta = set()
	for seq in all_seqs:
		species_name = getPropertyFromFastaSeqHeader(seq.description, "organism")
		species_names_from_fasta.add(species_name)
	species_names_from_fasta = sorted(species_names_from_fasta)  # TODO: check how come names has no duplicates
	context.species_list_names = list(species_names_from_fasta)

	logger.debug("species_names_from_fasta")
	logger.debug(species_names_from_fasta)
	logger.debug("context.species_list_names")
	logger.debug(context.species_list_names)

	# Estimation of running time:
	with open(context.summary_file,'a') as f_sum:
		if os.path.exists(context.largeSpeciesFileList):
			num_lines = sum(1 for line in open(context.largeSpeciesFileList))
		else:
			num_lines=0
		spc_num=num_lines+len(species_names_from_fasta)
		f_sum.write("Species with data in Genbank: %d\n" %spc_num) # check where 'Species list contains'

	time_estimation_for_user(context)
	logger.debug("!! merge_syn_seqs_with_accepted_species : context.syn_acc_TaxID_dict")
	logger.debug(context.syn_acc_TaxID_dict)
	logger.debug("!! merge_subsp_seqs_with_parent_species:   context.parent_species_names_by_subsp_tax_id")
	logger.debug(context.parent_species_names_by_subsp_tax_id)
#
	#Merge Synonym names to their Accepted names:
	if context.UserFlags_dict['NameResType'] == 'MatchedNames' or context.UserFlags_dict['NameResType'] == 'SpellingCorrection':
		logger.debug("User selection: Merge synonyms to their accepted species, according to user's table") #in plants we can provide out table
		merge_syn_seqs_with_accepted_species(context)
	else:
		all_seqs = list(SeqIO.parse(context.fasta_seq_filename, "fasta"))
		already_in_list=[]
		for seq in all_seqs:
			original_tax_id = getPropertyFromFastaSeqHeader(seq.description, "taxonid")
			original_name = getPropertyFromFastaSeqHeader(seq.description, "organism")
			#Update original vs acc name dict:
			logger.info("Merge Syn was not selected, initialize array for Merge Subsp:")
			context.original_vs_newName_dict[str(original_name)]=str(original_name)
			context.spc_tax_afterMatch_dict[str(original_name)]=str(original_tax_id)

	for name in context.original_vs_newName_dict.values():
		if 'subsp.' in name or 'var.' in name:
			if 'subsp.' in name: cut_str = ' subsp.'
			else: cut_str = ' var.'
			name_binomial = name.split(cut_str)[0]
			logger.debug("Found subspecies taxa %s, check if binomial name %s in list" %(name,name_binomial))
			#Check if binomial name already in the dictionary so we wouldn't have replicas:
			if name_binomial in context.parent_species_names_by_parent_tax_id.values():
				parent_taxId=get_parent_taxId(name_binomial,context.parent_species_names_by_parent_tax_id)
				if parent_taxId == 'Failed':
					return 0
				context.parent_species_tax_id_by_subsp_tax_id[context.spc_tax_afterMatch_dict[name]]=parent_taxId
				context.parent_species_names_by_parent_tax_id[parent_taxId]=name_binomial
			#Check if binomial name in the list of spc so it will be merged:
			elif name_binomial in context.original_vs_newName_dict.values():
				logger.debug("Found binomial name, add to merge dictionaries")
				context.parent_species_tax_id_by_subsp_tax_id[context.spc_tax_afterMatch_dict[name]]=context.spc_tax_afterMatch_dict[name_binomial]
				context.parent_species_names_by_parent_tax_id[context.spc_tax_afterMatch_dict[name_binomial]]=name_binomial
			#Only replace the name
			else:
				logger.debug("Replace name only")
				context.parent_species_tax_id_by_subsp_tax_id[context.spc_tax_afterMatch_dict[name]]=context.spc_tax_afterMatch_dict[name]
				context.parent_species_names_by_parent_tax_id[context.spc_tax_afterMatch_dict[name]]=name_binomial
		else:
			context.parent_species_tax_id_by_subsp_tax_id[context.spc_tax_afterMatch_dict[name]]=context.spc_tax_afterMatch_dict[name]
			context.parent_species_names_by_parent_tax_id[context.spc_tax_afterMatch_dict[name]]=name

	if context.UserFlags_dict['Merge_Subsp'] == 'on':
		logger.debug("User selection: Merge subsp. amd variants to their parents species")
		merge_subsp_seqs_with_parent_species(context)

	logger.debug(len(context.cluster_contexts))

	# ------------------------------------

	if not verify_cluster_data_is_enough(context, init_check=True):
		statusFail_LogFile(context, 'Less than 5 species after merge operation')
		with open(context.working_dir+'/calc_time_vars.txt','w') as f_time_calc:
			f_time_calc.write(' ')
			f_time_calc.close()
		return

	# ------------------------------------

	cache_prefix = "clust_"
	dirs_for_cache = [context.clustering_dir, context.cluter_its_dir, context.clustering_results_dir]
	cluster_seqs(cluster_method, cluster_script_dir, id, context)
	cache_file(context.fasta_seq_filename, dirs_for_cache, id, cache_prefix)
	final_clusters_dir = context.clustering_results_dir

	all_fasta_file_list = [os.path.join(context.clustering_results_dir, f) for f in
						   os.listdir(context.clustering_results_dir)]
	all_fasta_file_list = sorted(all_fasta_file_list, key=lambda fasta: get_number_of_different_species(fasta),
								 reverse=True)

	for cluster_file in all_fasta_file_list:
		diff_taxa_in_cluster = get_number_of_different_species(cluster_file)
		if diff_taxa_in_cluster < int(ott_config['general']['min_species_in_cluster']):
			logger.info("Less than min number of species in: %s" % cluster_file)
		else:
			logger.info("Adding cluster for file %s" % cluster_file)
			context.add_cluster(cluster_fasta_file=cluster_file)

	# In case there is #
	if os.stat(os.path.join(context.cluter_its_dir, context.id + "-allseq-its-only_.fasta")).st_size == 0:
		logger.debug("No ITS sequences for this run")
	elif context.its_min_taxa == 'no':
		logger.debug("Not enough ITS sequences for this run")
	else:
		if context.its_support:
			logger.info("Adding ITS FASTA files as ITS clusters")
			its_context = context.add_its_final_cluster()
			#Create ITS flagFile:
			with open(its_context.cluter_work_dir + '/ITS_CLUSTER', 'w') as f_flag_its:
				f_flag_its.close()
			if its_context is None:
				logger.info("ITS final cluster was empty and was not added")

	# ------------------------------------

	if not verify_cluster_data_is_enough(context, init_check=True):
		statusFail_LogFile(context, 'Less than 5 species after merge operation')
		return

	# ------------------------------------


	context.optimize_active_clusters()

	# MD- Need to add -> Change the code so it will run in parallel:
	for cluster_context in context.cluster_contexts:
		cluster_dir = os.path.dirname(cluster_context.all_seqs_fasta_filename)
		if not os.path.exists(cluster_dir+'/ITS_CLUSTER'):
			init_cluster_mult_acc_and_db(cluster_context, context, cluster_script_dir, outgrouup_selection)
		else:
			logger.debug("This is an ITS cluster, skip blast and copy file")
			shutil.copyfile(cluster_context.all_seqs_fasta_filename,cluster_dir + '/seqs-no-multiple-accessions.fasta')

	# Sorting the cluster contexts according to the number of species in each cluster
	context.cluster_contexts.sort(key=lambda cluster_context: cluster_context.get_number_of_different_species(),
								  reverse=True)

	context.init_cluster_contexts_for_loci_trees_list()

	if not verify_cluster_data_is_enough(context):
		statusFail_LogFile(context, 'Not enough data in clusters')
		return


	#---------------------------------------------------------------------
	# NEW NEW NEW
	#MD (NEW)  -  Check if this is the correct position for adding Large Species: species with large number of sequences in genbank:
	#---------------------------------------------------------------------
	if os.path.exists(context.largeSpeciesFileList):
		Check_and_Add_LargeTaxIds(context)
	else:
		logger.debug("No Large species fasta file")

	#Copy to Summary Directory:
	if os.path.exists(context.largeSpeciesDir + '/LOG_blastall.txt'):
		shutil.copyfile(context.largeSpeciesDir + '/LOG_blastall.txt', context.summary_dir + '/LOG_blastall.txt')
		large_dir_summary = context.summary_dir + '/LargeSpecies'
		create_dir_if_not_exists(large_dir_summary)
		for cluster_num in Dict_Add_gi_to_Cluster.keys():
			cluster_to_add = str(cluster_num)
			gi_to_add = Dict_Add_gi_to_Cluster[cluster_num]
			shutil.copyfile(context.largeSpeciesDir + "add_seq_" + str(cluster_num) + '_.fasta',
							large_dir_summary + "add_seq_" + str(cluster_num) + '_.fasta')

	##############################################################################
	# Concatenating the alignments
	##############################################################################

	logger.info(
		"=======================================================================================================")
	logger.info(
		"====================== Perform concatenation of all clusters ====================================")
	logger.info(
		"=======================================================================================================")

	#Perform concat of all clusters, add outgroup sequence and align files:
	#process_seq_for_concat_tree(context, id,outgrouup_selection)
	process_seq_and_create_concat_file(context, id, outgrouup_selection)

	#Create Nexus MSA file:
	make_nexus_msa_and_web_partition(context)

	pickle_context(context, context.pickled_name, context.pickled_json_name)

	# Creating an empty file so that other processes can know that the job is done successfully
	logger.debug("Marking pretree process is done - creating %s " % context.pre_tree_done)
	open(context.pre_tree_done, 'a').close()

	#Added by Michal:
	# Verify concat file exists:
	if os.path.exists(context.concat_seqs_fasta_filename):
		if os.stat(context.concat_seqs_fasta_filename).st_size == 0:
			logger.info("Concat File exists but Empty: %s" % context.concat_seqs_fasta_filename)
			#Result indication file:
			f_end_status_name = context.summary_dir + '/AlignEnd_FinishedExecutionWithEmptyConcatFile'
			f_end = open(f_end_status_name, 'w')
			f_end.close()
			statusFail_LogFile(context, 'Alignment ended with an empty Concat file')
			return
		else:
			logger.info("Concat File exists: %s" % context.concat_seqs_fasta_filename)
			#Result indication file:
			#with open(context.final_status,'a') as f_status:
			#	f_status.write("Alignment Finished successfully with a concat file\n")
			f_end_status_name = context.debug_dir + '/AlignEnd_FinishedExecutionWithConcatFile'
			f_end = open(f_end_status_name, 'w')
			f_end.close()
	else:
		logger.info("Missing concat file: %s" % context.concat_seqs_fasta_filename)
		#Result indication file:
		f_end_status_name = context.summary_dir + '/AlignEnd_FinishedExecutionMissingConcatFile'
		f_end = open(f_end_status_name, 'w')
		f_end.close()
		statusFail_LogFile(context, 'Alignment ended but Concat file is Missing')
		return

	#Copy alignment file to Summary dir:
	shutil.copyfile(context.concat_seqs_fasta_filename,
					context.summary_dir + '/' + context.id + "-concat-aligned.fasta")

	#Check if user chose clusters data only:
	if str(context.OTTOptions) == 'clustering':
		logger.info(
			"=======================================================================================================")
		logger.info(
			"====                     Finished Clustering and Alignment execution                              =====")
		logger.info(
			"====                                   D O N E                                                    =====")
		logger.info(
			"=======================================================================================================")
		#Zip_summary dir:
		create_summary_zip(context)
		#f_fakeTree = open(context.summary_dir + '/Result_Tree_' + context.id + '.tre', 'w')
		#f_fakeTree.close()
		job_pass = open(context.working_dir + '/JOB_PASS.txt', 'w')
		job_pass.close()
		with open(context.final_status, 'a') as f_status:
			f_status.write("PASS - Clusters data created successfully\n")
		return

	logger.info(
		"=======================================================================================================")
	logger.info(
		"====                           Finished Alignment execution                                       =====")
	with open(context.debug_dir+'/msa_concat_done.txt','w') as f_concat_file_done:
		f_concat_file_done.close()
	logger.info(
		"====                              Start Phylogenetic tree                                         =====")


	#Build_MB_Tree(f_end_status_name,context)
	Create_Tree(context,f_end_status_name)


if __name__ == "__main__":

	start_time = time.time()

	print("-------------------------------------in the begining of the file----------------------------")
	parser = argparse.ArgumentParser(description='Main ploiDB script - builds phylogenetic trees per genus')
	parser.add_argument('--taxa-list-file', '-n', help='File with the list of taxa to run the pipeline on',
						required=True)
	parser.add_argument('--working-dir', '-w', help='File with the list of genera to run the pipeline on',
						required=False, default=None)
	#parser.add_argument('--log-filename', '-l', help='log filename', required=False, default=None)
	#parser.add_argument('--debug-filename', '-d', help='log filename', required=False, default=None)
	#parser.add_argument('--cluster-method', '-m', help='external scripts dir', required=True)
	parser.add_argument('--params-path', '-pp', help='Params file', required=False)
	parser.add_argument('--config-filename', '-c', help='Config file', required=True)
	parser.add_argument('--id', '-i', help='ID of run', required=True)

	parser.add_argument('--include-subsp', dest='ignore_subsp',
						help='boolean flag if should include subspecies and var.',
						required=False, action='store_false', default=False)
	parser.add_argument('--ignore-subsp', dest='ignore_subsp', help='boolean flag if should ignore subspecies and var.',
						required=False, action='store_true')
	parser.add_argument('--merge-subsp', dest='merge_subsp',
						help='boolean flag if should merge subspecies and var. into species',
						required=False, action='store_true', default=False)

	args = parser.parse_args()

	print("run arguments:")
	dict_args = vars(args)

	for key in dict_args:
		print("%s : %s" % (key, dict_args[key]))

	# reading all taxon names from taxon names file
	file = open(args.taxa_list_file, 'r')
	taxa_list = file.read().splitlines()
	taxa_list = [name.strip() for name in taxa_list]
	taxa_list = [name.replace('_',' ') for name in taxa_list]
	file.close()


	id = args.id
	working_dir = args.working_dir
	log_filename = working_dir + '/Log-' + id + '.txt'
	debug_filename = working_dir + '/Debug-' + id + '.txt'
	cluster_method = 'orthomclits'
	should_ignore_subsp = args.ignore_subsp
	should_merge_subsp = args.merge_subsp

	#Check if params file exists, otherwise use default file :
	#try:
	create_dir_if_not_exists(working_dir)
	if args.params_path:
		shutil.copyfile(args.params_path, working_dir + '/params.txt')
		logger.debug("User params file : %s" % args.params_path)

	#try:
	generate_taxa_tree(id, taxa_list, working_dir, log_filename, debug_filename, cluster_method, should_ignore_subsp,
						   should_merge_subsp)

	#except: # catch *all* exceptions
	#	create_dir_if_not_exists(working_dir + '/SummaryDir')
	#	logger.debug("Exception error on generate_taxa_tree!!!!")
	#	with open(working_dir + '/SummaryDir/FianlStatus.txt','w') as sum_Status_f:
	#		sum_Status_f.write("OTT ERROR. Please contact us with your job ID number")

	#running time calc:
	total_time_min = (time.time() - start_time) / 60

