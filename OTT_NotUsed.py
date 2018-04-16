import datetime
import time
from PloiDbContext import PloiDbContext
from ploidbCommon import *
from blastSequencesLocaly import *
#from getGenusSequences import get_taxa_sequences, get_species_ids_and_names_from_taxa_list
from getGenusSequences import *
from ClusterContext import SeqType
from Create_Tree import *		# Tree generation functions: part2_create_mrbayes_blk, create_mb_config_file
import logging.handlers
import mysql.connector
import os
import re
import json
import zipfile
import argparse
from ctypes import *
from ploidbUtils import *
from PathHelper import PathHelper
from name_resolve import do_resolve_names
from Bio import Phylo
from Bio import AlignIO
from handleMultAccessions import handle_multiple_accessions
from revClusteredFasta import reverse_seqs_in_fasta
from adjustDirFasta import adjust_direction_fasta_using_mafft
from MBconfig_Dictionary import mb_block_dict
import createJobFile
import pandas as pd
from ete3 import NCBITaxa

__author__ = 'Michal'


def create_locus_tree(cluster_context, context):
	if cluster_context.is_its_cluster:
		logger.info("Calling mafft in order to add ITS outgroup to the sequence fasta")
		add_seq_to_existing_alignment(cluster_context.all_seqs_fasta_filename_no_multiple_accessions,
									  cluster_context.outgroup_seq_full_filename,
									  cluster_context.all_seqs_with_outgroup_filename, cluster_context.cluter_work_dir)
		shutil.copyfile(cluster_context.all_seqs_with_outgroup_filename,
						cluster_context.all_seqs_aligned_fasta_filename)
		logger.info("DONE Adding ITS outgroup")
	else:
		merge_seq_files(
			[cluster_context.all_seqs_fasta_filename_no_multiple_accessions,
			 cluster_context.outgroup_seq_full_filename],
			cluster_context.all_seqs_with_outgroup_filename,'False')
		logger.info("DONE Writing outgroup")

		logger.info("START handling reversed seqs")
		adjust_direction_fasta_using_mafft([cluster_context.all_seqs_with_outgroup_filename])
		# reverse_seqs_in_fasta(cluster_context,genus_name)
		logger.info("DONE handling reversed seqs")

		logger.info("Aligning %s. output will be written to %s" % (
			cluster_context.all_seqs_with_outgroup_filename, cluster_context.all_seqs_aligned_fasta_filename))

		if os.path.exists(cluster_context.all_seqs_aligned_fasta_filename):
			logger.debug("File %s already exists - skipping GUIDANCE" % cluster_context.all_seqs_aligned_fasta_filename)
		elif not context.should_run_guidance:
			logger.info("run guidance option is off - skipping GUIDANCE")
			logger.debug("context.should_run_guidance is FALSE, copying file %s to %s" %
						 (cluster_context.all_seqs_with_outgroup_filename,
						  cluster_context.all_seqs_aligned_fasta_filename))
			shutil.copyfile(cluster_context.all_seqs_with_outgroup_filename,
							cluster_context.all_seqs_aligned_fasta_filename)
		else:
			logger.debug("File %s not found - calling GUIDANCE" % cluster_context.all_seqs_aligned_fasta_filename)
			run_guidance_for_seqs_and_columns(cluster_context.all_seqs_with_outgroup_filename,
											  cluster_context.guidance_work_dir1,
											  cluster_context.guidance_work_dir2,
											  cluster_context.all_seqs_aligned_fasta_filename, context)

	logger.info("Rewriting %s to %s - using organism name as seq description" % (
		cluster_context.all_seqs_aligned_fasta_filename, cluster_context.all_seqs_for_tree_creation_fasta_filename))
	names_dict1 = rewrite_fasta_with_organism_name_as_desc(cluster_context.all_seqs_aligned_fasta_filename,
														   cluster_context.all_seqs_for_tree_creation_fasta_filename)
	names_dict2 = rewrite_fasta_with_organism_name_as_desc(cluster_context.outgroup_seq_full_filename,
														   cluster_context.outgroup_seq_filename)

	logger.info("Saving name mapping in %s" % cluster_context.rewrite_names_dict_filename)
	all_names_dict = names_dict1.copy()
	all_names_dict.update(names_dict2)
	with open(cluster_context.rewrite_names_dict_filename, 'w') as names_handle:
		json.dump(all_names_dict, names_handle)

	logger.info("DONE Rewriting")

def init_common_NEED_TO_CHECK_outgroup(ploidb_context):
	# Open file with Outgroup outcome - either None or selected outgroup name:
	f_outgroup = open(ploidb_context.outgroup_file, 'w')
	# This way we ensure that the clusters table exists
	write_stats_to_db(ploidb_context)


	with get_db_conn(ploidb_context.db_file) as conn:
		curs = conn.cursor()

		# Iterating all possible outgroups - the WHERE condition is only optimization to make sure very small outgroups are not checked
		#query = "select outgroup_name,outgroup_type,num_of_clusters,data_matrix_size " \
		#		"from outgroup_stats " \
		#		"where data_matrix_size >= 0.25 * (select max(data_matrix_size) from outgroup_stats) limit 50"

		#query = "select outgroup_name,outgroup_type,num_of_clusters,data_matrix_size " \
		#		"from outgroup_stats " \
		#		"where data_matrix_size >= 0.25 * (select max(data_matrix_size) from outgroup_stats) limit 50 " \
		#		"AND outgroup_name NOT LIKE '%%s%'"
#
		#curs.execute(query)
		#outgroup_rows = curs.fetchall()
		#logger.info("Found %d candidates for outgroup" % len(outgroup_rows))
		#get all genera in current Alignment species list:
		genera_list=[]
		for specie in ploidb_context.species_list_names:
			specie_split=specie.split(' ')
			#Add only genus not already in list:
			if specie_split[0] not in genera_list:
				genera_list.append(specie_split[0])
			logger.debug("Init Common Outgroup, genera to exclude: %s\n" % specie_split[0])
		outgroup_context_list = list()
		get_next_outgroup=0
		query_str=''
		#get results that will not include names
		for genus in genera_list:
			#create query string:
			query_str+= ("AND outgroup_name NOT LIKE '%" + genus + "%' ")
			query = "select outgroup_name,outgroup_type,num_of_clusters,data_matrix_size " \
					"from outgroup_stats " \
					"where data_matrix_size >= 0.25 * (select max(data_matrix_size) from outgroup_stats) limit 50 %s" %query_str \


		curs.execute(query)
		outgroup_rows = curs.fetchall()
		logger.info("Found %d candidates for outgroup" % len(outgroup_rows))

		for outgroup_row in outgroup_rows:
			outgroup_name = outgroup_row[0]
			#Check if legal outgroup
			print("%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check outgroup %s" % outgroup_name)
			print("%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check outgroup %s" % outgroup_name)
			if (Check_outgroup_forLegalName(outgroup_name) == 1):
				logger.debug("Remove %s from outgroup list, illegal name (No Name resolution match)" % (outgroup_name))
				continue
			#Check if in Alignment list:
			for genus in genera_list:
				if genus in outgroup_name:
					logger.debug("Remove %s from outgroup list, includes genus name %s" % (outgroup_name,genus))
					get_next_outgroup=1
					break
			if get_next_outgroup == 1:
				get_next_outgroup=0
				continue
			else:
				get_next_outgroup=0
			outgroup_type = outgroup_row[1]
			num_of_clusters = outgroup_row[2]
			data_matrix_size = outgroup_row[3]

			logger.debug("@@@@@@@@ Processing the outgroup %s @@@@@@@" % outgroup_name)
			outgroup_context = get_best_clusters_coverage_for_outgroup(conn, ploidb_context, outgroup_name,
																	   outgroup_type, num_of_clusters, data_matrix_size)
			outgroup_context_list.append(outgroup_context)

		#Check if outgroup_list is Empty:
		if not outgroup_context_list:
			logger.info("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
			logger.info("^^^^^^^^^                      Outgroup list is EMPTY !!!                 ^^^^^^^^^^^^^^^")
			logger.info("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
			f_outgroup.write('None')
			f_outgroup.close()
			return False

		selected_outgroup_context, top_outgroup_contexts = get_best_outgroup_context(outgroup_context_list)
		logger.info("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
		logger.info(
			"^^^^^^^^^ Outgroup selected is best_outgroup_context=%s ^^^^^^^^^^^^^^^" % selected_outgroup_context)
		logger.info("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
		f_outgroup.write('%s' % selected_outgroup_context)
		f_outgroup.close()
		ploidb_context.cluster_contexts_for_concat_tree = ploidb_context.get_clusters_by_ids(
			selected_outgroup_context.get_selected_cids())
		ploidb_context.top_outgroup_contexts = top_outgroup_contexts
		ploidb_context.selected_outgroup = selected_outgroup_context

		for cluster_context in ploidb_context.cluster_contexts_for_concat_tree:
			cluster_context.is_used_for_concat_tree = True

		Write_Clusters_data_File(ploidb_context)
		write_outgroup_for_concat_to_file(ploidb_context, selected_outgroup_context)


def process_seq_for_concat_tree(context, genus_name,outgrouup_selection):

	if outgrouup_selection == 'None':
		logger.info("No Outgroup selection for  " + genus_name)
		init_NO_common_outgroup(context,outgrouup_selection)
	else:
		logger.info("Initializing common outgroup  for  " + genus_name)
		init_common_outgroup(context)
		# For Outgroup option, check if seqs are ready:
		if is_seqs_for_concat_tree_exists(context):
			logger.info("All sequences for concat tree creating exist. Skipping their creation")
		else:
			logger.info("Writing seq files that contains the concat outgroup")
			for idx, cluster_context in enumerate(context.cluster_contexts_for_concat_tree):
				if cluster_context.is_its_cluster:
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
					logger.info("Writing seq file that contains the concat outgroup for into %s" % cluster_context.all_seqs_with_concat_outgroup_filename)
					if cluster_context.no_of_outgroup_sequence_for_concat > 0:
						merge_seq_files([cluster_context.all_seqs_fasta_filename_no_multiple_accessions,
										 cluster_context.outgroup_sequence_for_concat_filename],
										cluster_context.all_seqs_with_concat_outgroup_filename,'False')
					else:
						shutil.copyfile(cluster_context.all_seqs_fasta_filename_no_multiple_accessions,
										cluster_context.all_seqs_with_concat_outgroup_filename)
					logger.info("START handling reversed seqs for concat")
					adjust_direction_fasta_using_mafft([cluster_context.all_seqs_with_concat_outgroup_filename])
					logger.info("DONE handling reversed seqs for concat")
					fasta_for_alignment = cluster_context.all_seqs_with_concat_outgroup_filename

					if os.path.exists(cluster_context.all_seqs_concat_aligned_fasta_filename):
						logger.debug("File %s already exists - skipping GUIDANCE" % cluster_context.all_seqs_concat_aligned_fasta_filename)
					else:
						logger.debug("3GUIDANCE Flag: %s" %context.should_run_guidance)
						if context.should_run_guidance != 'True':
							logger.info("run guidance option is off - skipping GUIDANCE")
							logger.debug("context.should_run_guidance is FALSE, copying file %s to %s" %
										 (cluster_context.all_seqs_with_concat_outgroup_filename,
										  cluster_context.all_seqs_concat_aligned_fasta_filename))
							run_alignment(fasta_for_alignment,cluster_context.all_seqs_concat_aligned_fasta_filename,context,cluster_context.cluster_concat_work_dir)
							#shutil.copyfile(cluster_context.all_seqs_with_concat_outgroup_filename,
							#				cluster_context.all_seqs_concat_aligned_fasta_filename)
						else:
							logger.debug(
								"File %s not found - calling GUIDANCE" % cluster_context.all_seqs_concat_aligned_fasta_filename)
							run_guidance_for_seqs_and_columns(fasta_for_alignment, cluster_context.guidance_concat_work_dir1,
															  cluster_context.guidance_concat_work_dir2,
															  cluster_context.all_seqs_concat_aligned_fasta_filename, context)

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

	if outgrouup_selection == 'None':
		concat_alignment(context)
	else: # outgrouup_selection == 'Single':
		concat_alignment_SingleOutgroup(context)

	logger.info("DONE Concatenating alignments for  " + genus_name)


def init_user_common_outgroup(ploidb_context):
	# Open file with Outgroup outcome - either None or selected outgroup name:
	f_outgroup = open(ploidb_context.outgroup_file, 'w')
	# This way we ensure that the clusters table exists
	write_stats_to_db(ploidb_context)

	outgroup_context_list = list()
	get_next_outgroup=0

	##User outgroup is used:
	#outgroup_name = ploidb_context.User_Outgroup
	#logger.debug("@@---------------------------------------------------------------@@")
	#logger.debug("@@@@@@@@ User Outgroup selection is: %s" %outgroup_name)
	#logger.debug("@@---------------------------------------------------------------@@")

	outgroup_context = OutGroupSelectionContext(ploidb_context, outgroup_name, outgroup_type,
												clusters_ids_with_outgroup_seq, clusters_ids_without_outgroup_seq)

	#logger.debug("@@@@@@@@ Processing the outgroup %s @@@@@@@" % outgroup_name)
	#outgroup_context = get_best_clusters_coverage_for_outgroup(conn, ploidb_context, outgroup_name,
	#														   outgroup_type, num_of_clusters, data_matrix_size)
	#outgroup_context_list.append(outgroup_context)

	#Check if outgroup_list is Empty:
	if not outgroup_context_list:
		logger.info("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
		logger.info("^^^^^^^^^                      Outgroup list is EMPTY !!!                 ^^^^^^^^^^^^^^^")
		logger.info("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
		f_outgroup.write('None')
		f_outgroup.close()
		return False

	selected_outgroup_context, top_outgroup_contexts = get_best_outgroup_context(outgroup_context_list)
	logger.info("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
	logger.info(
		"^^^^^^^^^ Outgroup selected is best_outgroup_context=%s ^^^^^^^^^^^^^^^" % selected_outgroup_context)
	logger.info("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
	f_outgroup.write('%s' % selected_outgroup_context)
	f_outgroup.close()
	ploidb_context.cluster_contexts_for_concat_tree = ploidb_context.get_clusters_by_ids(
		selected_outgroup_context.get_selected_cids())
	ploidb_context.top_outgroup_contexts = top_outgroup_contexts
	ploidb_context.outgroupSelection = get_selected_outgroup(ploidb_context)
	with open(ploidb_context.summary_file,'w') as f_sum:
		f_sum.write("Selected Outgroup: %s\n" % ploidb_context.outgroupSelection)

	for cluster_context in ploidb_context.cluster_contexts_for_concat_tree:
		cluster_context.is_used_for_concat_tree = True

	Write_Clusters_data_File(ploidb_context)
	write_outgroup_for_concat_to_file(ploidb_context, selected_outgroup_context)


def RERUN_process_seq_and_create_concat_file(context, genus_name,outgrouup_selection):

	if outgrouup_selection == 'None':
		logger.info("No Outgroup selection for  " + genus_name)
		init_NO_common_outgroup(context,outgrouup_selection)
	else:
		logger.info("Initializing common outgroup  for  " + genus_name)
		init_common_outgroup(context)

	# For Outgroup option, check if seqs are ready:
	original_working_dir = context.working_dir  #context.UserFlags_dict['OriginJobDir']
	original_concat_dir = original_working_dir + '/concat/'
	logger.info("Writing seq files that contains the concat outgroup")
	#for idx, cluster_context in enumerate(context.cluster_contexts_for_concat_tree):
	for idx in context.rerun_chosen_clusters:
		f_all_seqs_fasta_filename_no_multiple_accessions = original_working_dir + '/' + str(idx) + "/seqs-no-multiple-accessions.fasta"
		f_all_seqs_concat_aligned_fasta_filename = original_concat_dir + str(idx) + "/seqs-aligned-concat.fasta"
		f_outgroup_sequence_for_concat_filename = original_concat_dir + str(idx)  + "/outgroup_seq_concat.fasta"
		f_all_seqs_with_concat_outgroup_filename = original_concat_dir + str(idx)  + "/seqs-with-out-group-concat.fasta"
		dir_cluster_concat_work_dir = original_concat_dir + str(idx)
		f_all_seqs_for_concat_tree_creation_fasta_filename = original_concat_dir + str(idx) + "/seqs-organism-concat.fasta"
		f_outgroup_sequence_for_concat_short_desc_filename = original_concat_dir + str(idx) + "/outgroup_seq_concat_shortdesc.fasta"

		#ITS cluster only:
		test_its= False
		if test_its:
			if outgrouup_selection == 'None':
					shutil.copyfile(f_all_seqs_fasta_filename_no_multiple_accessions,
								f_all_seqs_concat_aligned_fasta_filename)
			else:
				#in case of an outgroup, add it using mafft:
				if os.path.exists(f_all_seqs_with_concat_outgroup_filename):
					logger.info("Calling mafft in order to add ITS outgroup to the concat sequence fasta")
					add_seq_to_existing_alignment(f_all_seqs_fasta_filename_no_multiple_accessions,
											  f_outgroup_sequence_for_concat_filename,
											  f_all_seqs_with_concat_outgroup_filename,
											  dir_cluster_concat_work_dir)
					shutil.copyfile(f_all_seqs_with_concat_outgroup_filename,
								f_all_seqs_concat_aligned_fasta_filename)
				else:
					logger.info("No outgroup - copying %s to %s" % (
						f_all_seqs_fasta_filename_no_multiple_accessions,
						f_all_seqs_concat_aligned_fasta_filename))
					shutil.copyfile(f_all_seqs_fasta_filename_no_multiple_accessions,
									f_all_seqs_concat_aligned_fasta_filename)
				logger.info("DONE Adding ITS outgroup")
		else:
			#Handling other clusters, part of ITS:
			logger.info("Writing seq file that contains the concat outgroup for into %s" % f_all_seqs_with_concat_outgroup_filename)
			if outgrouup_selection == 'None':
				shutil.copyfile(f_all_seqs_fasta_filename_no_multiple_accessions,
									f_all_seqs_with_concat_outgroup_filename)
			else:
				# Clusters with outgroup:
				#if os.path.exists(f_all_seqs_with_concat_outgroup_filename):
				if os.path.exists(f_outgroup_sequence_for_concat_filename):
					merge_seq_files([f_all_seqs_fasta_filename_no_multiple_accessions,
									 f_outgroup_sequence_for_concat_filename],
									f_all_seqs_with_concat_outgroup_filename,'False')
				else:
					shutil.copyfile(f_all_seqs_fasta_filename_no_multiple_accessions,
									f_all_seqs_with_concat_outgroup_filename)
			logger.info("START handling reversed seqs for concat")
			adjust_direction_fasta_using_mafft([f_all_seqs_with_concat_outgroup_filename])
			logger.info("DONE handling reversed seqs for concat")
			fasta_for_alignment = f_all_seqs_with_concat_outgroup_filename

			#if os.path.exists(f_all_seqs_concat_aligned_fasta_filename):
			#	logger.debug("File %s already exists - skipping Alignment" % f_all_seqs_concat_aligned_fasta_filename)
			#else:
			# Perform Alignment:
			if context.UserFlags_dict['FilterMSA_Method'] == 'GUIDANCE':
				run_guidance_for_seqs_and_columns(fasta_for_alignment, dir_cluster_concat_work_dir+'/guidance1/',
													  dir_cluster_concat_work_dir+'/guidance2/',
													  f_all_seqs_concat_aligned_fasta_filename, context)
			else:
				run_alignment(fasta_for_alignment,f_all_seqs_concat_aligned_fasta_filename,context,dir_cluster_concat_work_dir)
				if context.UserFlags_dict['FilterMSA_Method'] is not 'None':
					# Check MSA Filter Methods:
					logger.debug("User selected Filter MSA Method is: %s" %context.UserFlags_dict['FilterMSA_Method'])
					perform_filter_msa(context,f_all_seqs_concat_aligned_fasta_filename)

		logger.info("Rewriting %s to %s - using organism name as seq description" % (
			f_all_seqs_concat_aligned_fasta_filename,
			f_all_seqs_for_concat_tree_creation_fasta_filename))
		names_dict1 = rewrite_fasta_with_organism_name_as_desc(
			f_all_seqs_concat_aligned_fasta_filename,
			f_all_seqs_for_concat_tree_creation_fasta_filename, False)
		names_dict2 = dict()
		#if cluster_context.no_of_outgroup_sequence_for_concat > 0:
		if os.path.exists(f_all_seqs_with_concat_outgroup_filename):
			names_dict2 = rewrite_fasta_with_organism_name_as_desc(
				f_outgroup_sequence_for_concat_filename,
				f_outgroup_sequence_for_concat_short_desc_filename, False)
		else:
			logger.debug("No outgroup found for cluster %s" % cluster_context)

		all_names_dict = names_dict1.copy()
		all_names_dict.update(names_dict2)
		with open(context.rewrite_names_dict_filename, 'w') as names_handle:
			json.dump(all_names_dict, names_handle)

		logger.info("DONE Rewriting")

	logger.info("DONE Writing seq files that contains the concat outgroup")

	#context.init_chloroplast_cluster()
	#context.init_mitochondrial_cluster()
	#context.init_nuc_cluster()

	logger.info("Concatenating alignments for  " + genus_name)
	#concat_alignment_SingleOutgroup(context)
	RERUN_concat_alignment_SingleOutgroup(context)

	logger.debug("context.outgroupSelection-----------------------------------------------------------------\n\n")
	logger.debug(context.outgroupSelection)

	logger.info("DONE Concatenating alignments for  " + genus_name)
	#Create_species_list included in the Alignment:
	with open(context.concat_seqs_fasta_filename,'r') as f_alignment:
		f_finalNames = open(context.final_species_file, 'w')
		for line in f_alignment:
			if (context.outgroupSelection) and ('>' in line) and (context.outgroupSelection not in line):
				line = line.strip()
				spc_name = line.replace('>','').replace('_',' ')
				f_finalNames.write(spc_name+'\n')
		f_finalNames.close()
		f_alignment.close()


def check_name_inNCBIdb(taxId):

	db_name = ploidb_config['general']['NCBI_NAMES_TAXIDS_DB']
	conn_ncbi_name = sqlite3.connect(db_name)
	curs_ncbi = conn_ncbi_name.cursor()

	query_name = "select Name from names_vs_tax_ids where TaxID is '%s'" % taxId
	logger.info("select Name from names_vs_tax_ids where TaxID is '%s'" % taxId)

	curs_ncbi.execute(query_name)
	rows_ncbi = curs_ncbi.fetchall()
	for line in rows_ncbi:	# query is empty
		#Handle sp. and other markers in names:
		spDot_remove = line[0].replace('sp.','sp')
		Bold_remove = spDot_remove.replace('BOLD:','BOLD-')
		coded_name = Bold_remove.replace('_',' ')
		return coded_name

##  # Ploidb regular version:
##	db_name = ploidb_config['general']['NCBI_TAXID_DB']
##	conn_ncbi_name = sqlite3.connect(db_name)
##	curs_ncbi = conn_ncbi_name.cursor()
##
##	query_name = "select Coded_Name from NCBI_AfterNameRes where Id is '%s'" % taxId
##	logger.info("select Coded_Name from NCBI_AfterNameRes where Id is '%s'" % taxId)
##
##	curs_ncbi.execute(query_name)
##	rows_ncbi = curs_ncbi.fetchall()
##	for line in rows_ncbi:	# query is empty
##		coded_name = line[0].replace('_',' ')
##		return coded_name


# TODO: better handle matching between names and ids
def resolve_names(context):

	ids_to_resolved_names = dict()
	ids_not_resolved = list()
	taxId_organism_dict=dict()

	logger.info("Resolving names")
	all_seqs = list(SeqIO.parse(context.fasta_seq_filename, "fasta"))

	organism_names = set()
	taxonids = set()
	for seq in all_seqs:
		organism_name = getPropertyFromFastaSeqHeader(seq.description, "organism")
		taxonid = getPropertyFromFastaSeqHeader(seq.description, "taxonid")
		organism_names.add(organism_name)
		taxonids.add(taxonid)
		taxId_organism_dict[organism_name]=str(taxonid)
	organism_names = sorted(organism_names)  # TODO: check how come names has no duplicates
	taxonids = sorted(taxonids)

	original_name_to_id = dict()
	with open(context.names_to_resolve, "w") as names_handle:
		writer = csv.writer(names_handle, delimiter=',', dialect=csv.excel)
		writer.writerow(["Id", "Name"])
		for index, organism_name in enumerate(organism_names):
			original_name_to_id[organism_name] = str(index)
			writer.writerow([index, organism_name])
			# Check if Accepted_name:
			id = str(index)
			# take Name from NCBI db and use this name for the search:
			name_organism_ncbi = check_name_inNCBIdb(taxId_organism_dict[organism_name])
			#if organism_name in context.names_vs_accepted_dict:

			#context.species_names_by_tax_id[tax_id] = accepted_noUnderscore
			checkTaxId = taxId_organism_dict[organism_name]
			if checkTaxId in context.species_names_by_tax_id.keys():
				newname = "%s" % context.species_names_by_tax_id[checkTaxId]
				ids_to_resolved_names[id] = newname
				logger.debug("Adding resolved name: %s ==> %s, organism %s, TaxId %s" % (id, newname, organism_name, checkTaxId))
#			if name_organism_ncbi in context.names_vs_accepted_dict:
#				newname = "%s" % context.names_vs_accepted_dict[name_organism_ncbi]
#				ids_to_resolved_names[id] = newname
#				logger.debug("Adding resolved name: %s ==> %s, organism %s" % (id, newname, name_organism_ncbi))
			else:
				ids_not_resolved.append(id)
				logger.debug("Removed organism, no accepted name according to genus: %s" % (name_organism_ncbi))
				pass  # TODO - remove the sequences from the list

	#sys.exit(0)
	resolved_seqs = list()
	resolved_names_taxonid = dict()  # since we need that the taxonid
	unresolved_names = set()
	for seq in all_seqs:
		# logger.debug("Processing seq.id=%s seq.description=%s" % (seq.id, seq.description))
		original_name = getPropertyFromFastaSeqHeader(seq.description, "organism")

		skip_seq = False
		resolved_organism_name = original_name
		if original_name in original_name_to_id:
			id = original_name_to_id[original_name]
			if id in ids_to_resolved_names and id not in ids_not_resolved:
				resolved_organism_name = ids_to_resolved_names[id]
			else:
				skip_seq = True
				if id not in ids_to_resolved_names:
					logger.info(
						"Filtering sequence due to name resolution - ID %s of %s was not found in dict for name resolution" % (
							id, original_name))
				if id in ids_not_resolved:
					logger.info(
						"Filtering sequence due to name resolution - ID %s of %s has a low or no name resolve score and will be removed" % (
							id, original_name))
		else:
			logger.error("ID was not found for Organism name %s" % original_name)

		if skip_seq:
			logger.info("Skipping seq of %s" % original_name)
			unresolved_names.add(original_name)
		else:
			original_taxonid = getPropertyFromFastaSeqHeader(seq.description, "taxonid")
			if resolved_organism_name not in resolved_names_taxonid:
				resolved_names_taxonid[resolved_organism_name] = original_taxonid

			#logger.debug("Renaming: %s ==> %s" % (original_name,resolved_organism_name))
			context.names_to_resolved_names[original_name] = resolved_organism_name
			original_genus = context.id
			original_genus = original_genus.lower().strip()
			resolved_genus = resolved_organism_name.split(" ")[0]
			resolved_genus = resolved_genus.lower().strip()

			logger.info("restandarize_seq, all_seq organism name resolve:  %s,%s to %s,%s" % (original_taxonid,original_name,resolved_names_taxonid[resolved_organism_name],resolved_organism_name))
			tax_id = resolved_names_taxonid[resolved_organism_name]
			restandarize_seq(seq, original_organism=original_name, organism=resolved_organism_name,
							 taxonid=resolved_names_taxonid[resolved_organism_name],
							 original_taxonid=original_taxonid)  # setting seq properties
			resolved_seqs.append(seq)

	shutil.copyfile(context.fasta_seq_filename, context.fasta_seq_org_names_filename)
	write_standarize_fasta(context.fasta_seq_filename, resolved_seqs)
	#sys.exit(0)
	context.unresolved_species_names = unresolved_names



#Get names and taxIds of species according to the taxId given:
#For Species
def get_species_name_and_ID(context, taxon_name,taxon_id, first_flag,should_ignore_subsp=False, should_merge_subsp=False):

	with open(context.species_names_taxid_file,mode="a",encoding='utf8',newline='') as names_taxid_f:
		writer = csv.writer(names_taxid_f,delimiter=',')
		if first_flag == 0: writer.writerow(['TaxID', 'Species_Name'])

		species_names_by_tax_id = dict()
		parent_species_names_by_subsp_tax_id = dict()
		parent_species_tax_id_by_subsp_tax_id = dict()
		writer.writerow([str(taxon_id), str(taxon_name)])
		context.species_names_by_tax_id[str(taxon_id)] = taxon_name
		names_taxid_f.close()

	names_taxid_f.close()
	return list(context.species_names_by_tax_id.keys()), list(context.species_names_by_tax_id.values()), parent_species_tax_id_by_subsp_tax_id, parent_species_names_by_subsp_tax_id
#For Genus
def get_Genus_species_name_and_ID(context, genus_name,first_flag,should_ignore_subsp=False, should_merge_subsp=False):

	with open(context.species_names_taxid_file,mode="a",encoding='utf8',newline='') as names_taxid_f:
		writer = csv.writer(names_taxid_f,delimiter=',')
		if first_flag == 0: writer.writerow(['TaxID', 'Species_Name'])

		genera_flag=0
		species_names_by_tax_id = dict()
		parent_species_names_by_subsp_tax_id = dict()
		parent_species_tax_id_by_subsp_tax_id = dict()

		# TPL db name:
		db_TPL_name = ploidb_config['general']['TPL_ALL_DB']
		conn_tpl = sqlite3.connect(db_TPL_name)
		curs_tpl = conn_tpl.cursor()
		query_species = "SELECT Specie FROM TPL_Name_vs_Accepted WHERE Accepted_Genus LIKE '%s' AND Status is not 'Unresolved'" % genus_name
		curs_tpl.execute(query_species)
		SpeciesUnderGenus_rows = curs_tpl.fetchall()
		specie_list_query=[]
		for line in SpeciesUnderGenus_rows:
			name = line[0].strip()
			specie_list_query.append(name)

		logger.info(specie_list_query)
		curs_tpl.close()

		# NCBI db name:
		db_NCBI_name = ploidb_config['general']['NCBI_NAMES_TAXIDS_DB']  # Names of Sientific Names and TaxIds (from names.dmp)
		conn_ncbi = sqlite3.connect(db_NCBI_name)
		curs_ncbi = conn_ncbi.cursor()
		species_list_forDB = Create_names_forDBselect(specie_list_query)
		query_taxId = "SELECT TaxID,Name FROM names_vs_tax_ids WHERE Name IN (%s)" % species_list_forDB
		logger.info(query_taxId)
		curs_ncbi.execute(query_taxId)
		SpeciesTaxIDs_UnderGenus_rows = curs_ncbi.fetchall()
		logger.info(SpeciesTaxIDs_UnderGenus_rows)
		logger.debug("Number of species found in data base for Genus %s: %d" % (genus_name,len(SpeciesTaxIDs_UnderGenus_rows)))
		if SpeciesTaxIDs_UnderGenus_rows:
			for line in SpeciesTaxIDs_UnderGenus_rows:
				writer.writerow(line)
				tax_id = str(line[0])
				species_name = str(line[1])
				context.species_names_by_tax_id[tax_id] = species_name
		else:
			logger.debug("Empty results for %s" %genus_name)
			return list(), list(), dict(), dict()

		curs_ncbi.close()

	names_taxid_f.close()

	return list(context.species_names_by_tax_id.keys()), list(context.species_names_by_tax_id.values()), parent_species_tax_id_by_subsp_tax_id, parent_species_names_by_subsp_tax_id

#For Family
def get_Family_species_name_and_ID(context, family_name,first_flag,should_ignore_subsp=False, should_merge_subsp=False):

	with open(context.species_names_taxid_file,mode="a",encoding='utf8',newline='') as names_taxid_f:
		writer = csv.writer(names_taxid_f,delimiter=',')
		if first_flag == 0: writer.writerow(['TaxID', 'Species_Name'])

		genera_flag=0
		species_names_by_tax_id = dict()
		parent_species_names_by_subsp_tax_id = dict()
		parent_species_tax_id_by_subsp_tax_id = dict()

		# TPL db name:
		db_TPL_name = ploidb_config['general']['TPL_ALL_DB']
		conn_tpl = sqlite3.connect(db_TPL_name)
		curs_tpl = conn_tpl.cursor()
		query_genus = "SELECT Genus FROM TPL_all WHERE Family LIKE '%s'" % family_name
		curs_tpl.execute(query_genus)
		GenusUnderFamily_rows = curs_tpl.fetchall()
		genera_list_query=[]
		for line in GenusUnderFamily_rows:
			genus_name = line[0].strip()
			genera_list_query.append(genus_name)

			query_species = "SELECT Specie FROM TPL_Name_vs_Accepted WHERE Accepted_Genus LIKE '%s' AND Status is not 'Unresolved'" % genus_name
			curs_tpl.execute(query_species)
			SpeciesUnderGenus_rows = curs_tpl.fetchall()
			specie_list_query=[]
			for line in SpeciesUnderGenus_rows:
				name = line[0].strip()
				specie_list_query.append(name)

			logger.info(specie_list_query)
			curs_tpl.close()

			# NCBI db name:
			db_NCBI_name = ploidb_config['general']['NCBI_NAMES_TAXIDS_DB']  # Names of Sientific Names and TaxIds (from names.dmp)
			conn_ncbi = sqlite3.connect(db_NCBI_name)
			curs_ncbi = conn_ncbi.cursor()
			species_list_forDB = Create_names_forDBselect(specie_list_query)
			query_taxId = "SELECT TaxID,Name FROM names_vs_tax_ids WHERE Name IN (%s)" % species_list_forDB
			logger.info(query_taxId)
			curs_ncbi.execute(query_taxId)
			SpeciesTaxIDs_UnderGenus_rows = curs_ncbi.fetchall()
			logger.info(SpeciesTaxIDs_UnderGenus_rows)
			logger.debug("Number of species found in data base for Genus %s: %d" % (genus_name,len(SpeciesTaxIDs_UnderGenus_rows)))
			if SpeciesTaxIDs_UnderGenus_rows:
				for line in SpeciesTaxIDs_UnderGenus_rows:
					writer.writerow(line)
					tax_id = str(line[0])
					species_name = str(line[1])
					context.species_names_by_tax_id[tax_id] = species_name
			else:
				logger.debug("Empty results for %s" %genus_name)
				#return list(), list(), dict(), dict()

			curs_ncbi.close()

	names_taxid_f.close()

	return list(context.species_names_by_tax_id.keys()), list(context.species_names_by_tax_id.values()), parent_species_tax_id_by_subsp_tax_id, parent_species_names_by_subsp_tax_id

# Converting taxon name to list of species (ids and names) using Entrez
def get_species_ids_and_names_from_taxon_name_NOT_IN_USE(context, taxon_name,first_flag,should_ignore_subsp=False, should_merge_subsp=False):

	names_taxid_f = open(context.species_names_taxid_file,mode="a",encoding='utf8',newline='')
	writer = csv.writer(names_taxid_f,delimiter=',')
	#Name,Authority,Original name,Original authority,Coded Name,Coded Authority,Score,Matched Name,Id
	#writer.writerow(['Name','Authority','Original name','Original authority','Coded Name','Coded Authority','Score','Matched Name','Id'])
	if first_flag==0: writer.writerow(['TaxID', 'Species_Name'])


	genera_flag=0
	species_names_by_tax_id = dict()
	parent_species_names_by_subsp_tax_id = dict()
	parent_species_tax_id_by_subsp_tax_id = dict()
	# TPL db name:
	db_name = ploidb_config['general']['NCBI_NAMES_TAXIDS_DB']  # Names of Sientific Names and TaxIds (from names.dmp)
	conn_tpl = sqlite3.connect(db_name)
	curs_tpl = conn_tpl.cursor()

	query_line,query_name_syn, genera_flag=CreateTPLdb_query(taxon_name,context)

	curs_tpl.execute(query_line)
	rows_from_db = curs_tpl.fetchall()
	# Check if there's a syn query:
	No_syn=0
	if query_name_syn != 'None':
		No_syn=1
		curs_tpl.execute(query_name_syn)
		names_syn_fromdb = curs_tpl.fetchall()


	if('TRUE'):
		for line in rows_from_db:
			tax_id = str(line[0])			# TaxId of the original name (not of the accepted name)
			species_name = str(line[1])		# Accepted name - According to Name resolution (should match Coded name without the '_'
			#verify no cf. or sp. in species names:
			if 'cf.' in species_name or 'sp.' in species_name:
				logger.debug("Skip species %s (include either cf. or sp.)" % species_name)
			else:
				writer.writerow(line)
				context.species_names_by_tax_id[tax_id] = species_name
			##original_name = str(line[2])	# Input for Name resolution (accepted and synonyms)
			##coded_name = str(line[4])		# Accepted Name with '_' -> should not contain any special signs
			##score = float(line[6])			# Score of Name resolution reliability
			##matched_name = str(line[7])		# Matched name in plant list - > synonym names as in plant list and Accepted as well (full names including authority): mainly for Tax Id recognition
			##tax_id = str(line[8])			# TaxId of the original name (not of the accepted name)
			##if score > 0.8 :			# Add printout -> Low score in name resolution, removed species XXX from list
			##	#species_dict[tax_id] = species_name
			##	writer.writerow(line)
			##	#writer.writerow([species_name,tax_id,original_name])
			##	context.CodedNames_TaxIds_dict[tax_id]=coded_name
			##	context.MatchedNames_TaxIds_dict[tax_id]=matched_name
			##	# NCBI name vs. Accepted species (Coded name):
			##	accepted_noUnderscore = coded_name.replace('_',' ')
			##	#context.names_vs_accepted_dict[species_name] = accepted_noUnderscore
			##	# Create TaxId match between Accepted and syn names:
			##	context.species_names_by_tax_id[tax_id] = accepted_noUnderscore
			##else:
			##	logger.debug('Removed Species: %s, TaxID %s (Name resolution score Lower than 0.8)' % (species_name,tax_id))
		#Check number of species:
		max_number_of_species = 20000
		number_of_species=len(context.species_names_by_tax_id)
		if number_of_species > max_number_of_species:
			logger.error("too many results for " + taxon_name + "(" + str(number_of_species) + "). Skipping")
			return list(), list()
		else:
			logger.debug("Number of species found in data base for Genus %s: %d" % (taxon_name,number_of_species))
			if number_of_species == 0:
				return list(), list(), dict(), dict()

		#Update Names and syn/name res file:
		if No_syn !=0:
			for line in names_syn_fromdb:
				name = str(line[0])
				accepted_name = str(line[1])
				#writer_syn.writerow([name,accepted_name])
				#context.names_vs_accepted_dict[name] = accepted_name


	names_taxid_f.close()
	#names_syn_f.close()

	#sys.exit(0)
	return list(context.species_names_by_tax_id.keys()), list(context.species_names_by_tax_id.values()), parent_species_tax_id_by_subsp_tax_id, parent_species_names_by_subsp_tax_id


#_____________________________________________________________________________________________________________________#

def get_children_of_parentTaxId(context, taxon_id_list):

	#NCBI Nodes db - connect abd get all kids
	children_taxId_list=[]
	potential_parents_list=[]
	nodes_DB = ploidb_config['general']['NCBI_Nodes_DB']  # Names of Sientific Names and TaxIds (from names.dmp)
	conn_ncbi = sqlite3.connect(nodes_DB)
	curs_ncbi = conn_ncbi.cursor()
	temp_taxa_list=[]
	for taxId in taxon_id_list:
		taxid_added = "'"+taxId+"'"
		temp_taxa_list.append(taxid_added)
	taxId_list_str = ','.join(temp_taxa_list)
	query_taxId = "SELECT tax_id FROM Nodes_dmp WHERE parent_tax_id IN (%s)" % taxId_list_str
	logger.debug(query_taxId)
	curs_ncbi.execute(query_taxId)
	all_children_lines = curs_ncbi.fetchall()
	temp_taxa_list=[]
	for line in all_children_lines:
		children_taxId_list.append(str(line[0]))
		new_taxId = "'"+str(line[0])+"'"
		temp_taxa_list.append(new_taxId)
	taxId_list_str = ','.join(temp_taxa_list)
	# get only taxIds that are parents, so we need to get their kids as well:
	query_taxId_whoRparents = "SELECT DISTINCT parent_tax_id FROM Nodes_dmp WHERE parent_tax_id IN (%s)"%taxId_list_str
	logger.debug(query_taxId_whoRparents)
	curs_ncbi.execute(query_taxId_whoRparents)
	potential_parents_lines = curs_ncbi.fetchall()
	for line in potential_parents_lines:
		potential_parents_list.append(str(line[0]))
	print(children_taxId_list)
	print(potential_parents_list)
	print(list(set(children_taxId_list) - set(potential_parents_list)))
	print(len(children_taxId_list))
	print(len(potential_parents_list))

	sys.exit()

	return children_taxId_list,potential_parents_list


# Michal: Get species taxids and names according to NCBI nodes db:
def get_species_TaxIds_and_Names(context,taxa_list, should_ignore_subsp=False, should_merge_subsp=False):
	logger.info("Converting taxa list to list of species ids and species names (%d)" %len(taxa_list))
	species_ids = []
	species_names = []
	parent_species_tax_id_by_subsp_tax_id = dict()
	parent_species_names_by_subsp_tax_id = dict()

	ncbi = NCBITaxa()
	#Check if digit -> TaxID taxon list:
	if taxa_list[0].isdigit():
		logger.info("Input taxa list was identified as TaxID list")
		switch_toTaxId = 0
	else:
		logger.info("Input taxa list was identified as Names list")
		switch_toTaxId = 1
	first_flag=0
	for taxon_item in taxa_list:
		taxon_item = taxon_item.rstrip()
		if switch_toTaxId == 1:
			nametranslator = ncbi.get_name_translator([taxon_item])
			if nametranslator:
				taxon_id = ncbi.get_name_translator([taxon_item])[taxon_item][0]
				logger.debug(taxon_id)
				logger.debug(ncbi.get_name_translator([taxon_item]))
			else:
				logger.debug("No TaxId for taxon %s" %taxon_item)
				continue
		else:
			taxon_id = taxon_item
		#Get all kids for taxID:
		#kids_taxIds_list, potential_parents_list = get_children_of_parentTaxId(taxon_id)    #taxId vs rank
		#remove taxIds who have no childeren:
		#remove_TaxIds_withNoKids(kids_taxIds_Dict)

	taxon_id_list=taxa_list
	get_children_of_parentTaxId(context,taxon_id_list)
	return species_ids, species_names, parent_species_tax_id_by_subsp_tax_id, parent_species_names_by_subsp_tax_id


def get_species_of_clade(clade_id):
	ncbi = NCBITaxa()
	species_list = ncbi.get_descendant_taxa(clade_id, collapse_subspecies=False,return_tree=True)
	#ncbi.get_descendant_taxa(4107, collapse_subspecies=False,return_tree=False)


def writeTaxonSequences_CdHit(context, fasta_filename, gb_filename, nextTaxId):
	temp_fasta_file = open(context.fasta_cdhit_temp_input_f, 'w')
	temp_fasta_o_file = open(context.fasta_cdhit_temp_output_f, 'w')
	temp_fasta_file.close()
	temp_fasta_o_file.close()

	# Check if hugeTree flag:
	hugeTree_flag = ploidb_config['general']['Huge_trees']

	#Check all paths:
#	for grp_name in GROUP_LIST:
#		sequences_path = getSpeciesSequenceFileLocation(nextTaxId,grp_name)  # getting gene data from the files downloaded from genebank
#		if sequences_path is None:
#			continue
#		else:
#			break

	sequences_path = getSpeciesSequenceFileLocation(nextTaxId)  # getting gene data from the files downloaded from genebank
	if sequences_path is None:
		# TODO: add code that checks if this is legitimate (i.e by querying NCBI online) or not
		# This can happen in cases where the taxon itself doesn't have sequences but its subspecies do
		logger.warning("Path not found for taxon " + str(nextTaxId) + ". Skipping this taxon")
		return 0
	else:
		#logger.info("getting sequences for " + nextTaxId + " from " + str(sequences_path))
		speciesSeqRecords = SeqIO.parse(sequences_path, "genbank")
		# The data is sorted so that on multiple runs of the pipeline, it would be easier to use cached results (which are based on file CRC)
		speciesSeqList = get_ordered_seqs(list(speciesSeqRecords))
		no_of_seq_for_taxonid = len(speciesSeqList)
		logger.debug("found " + str(no_of_seq_for_taxonid) + " sequences for " + nextTaxId)

		mac_seq_length = int(ploidb_config['general']['max_seq_length'])
		speciesSeqListWrittenToFasta = list()
		temp_fasta_file = open(context.fasta_cdhit_temp_input_f, 'w')
		for seq in speciesSeqList:
			if hugeTree_flag and no_of_seq_for_taxonid >= 1000:
				logger.info("HugeTree taxid %s exceeds Max number of seqs 1000" % (nextTaxId))
				no_of_seq_for_taxonid=0
				break
			if len(seq.seq) > mac_seq_length:
					logger.info("Seq length is %d for seq %s - ignoring this seq" % (len(seq.seq), seq.description))
			else:
				writeSequenceInFastaFormat(temp_fasta_file, nextTaxId, seq)
				speciesSeqListWrittenToFasta.append(seq)
		temp_fasta_file.close()
		if no_of_seq_for_taxonid > 1:
			logger.info('Running CD-HIT on Taxid - %s' % nextTaxId)
			CDHITcommand_to_exec = 'cd-hit -i %s -d 0 -o %s -c 0.98 -n 5 -G 1 -g 1 -b 25' % (context.fasta_cdhit_temp_input_f,context.fasta_cdhit_temp_output_f)
			exec_external_command_redirect_output(command_to_exec=CDHITcommand_to_exec,outfile=context.working_dir + "/CD-HIT_perTaxID.out",errfile=context.working_dir + "/CD-HIT_perTaxID.err")


		if fasta_filename is not None:
			with open(fasta_filename, 'a') as f1:
				if no_of_seq_for_taxonid > 1:
					read_cd_File = open(context.fasta_cdhit_temp_output_f, 'r')
				else:
					read_cd_File = open(context.fasta_cdhit_temp_input_f, 'r')
				lines = read_cd_File.readlines()
				for line in lines:
					if line.strip(): f1.write(line) # write only non-blank lines
			num_of_seq_AfterCdHit = (len(lines)/2)
			logger.debug("Wrote " + str(num_of_seq_AfterCdHit) + ' out of ' + str(no_of_seq_for_taxonid) + " new sequences in fasta format (after CD-Hit)")

		if gb_filename is not None:
			# logger.debug("Writing to GB file")
			all_genus_seq_gb_handle = open(gb_filename, "a")
			# TODO: check the sequece length and filter if its too long
			SeqIO.write(speciesSeqListWrittenToFasta, all_genus_seq_gb_handle, "genbank")
			all_genus_seq_gb_handle.close()
		# logger.debug("Wrote " + str(no_of_seq_for_taxonid) + " new sequences in genbank format")

		temp_fasta_file.close()
		temp_fasta_o_file.close()

		return no_of_seq_for_taxonid


# return_species_list_debug -> Working version 1 - just in case
def return_species_list_debug_wv1(input_list,context):
	# flags:
		# SpeciesDescendants:   if 'No' then we bring users species only and filter higher ranked list according to Filter flags (exact species list by user)
		#                       if 'Yes' then we bring also descendants of user input species and filter all according to Filter flags
	## or All_Descendants with exclude options:
		# Remove Hybrids
		# Remove Intraspecific (subsp, var...)
		# Remove Nomenclature (include: cf. / sp. / aff. / f.))
	# for Accepted names match need to prepare the correct dictionaries:
	# 1. for subsp merge:
	#			parent_tax_id = context.parent_species_tax_id_by_subsp_tax_id[original_tax_id]
	#			parent_name = context.parent_species_names_by_subsp_tax_id[original_tax_id]
	# 2. for Acc/Synonym merge:
	#			context.syn_acc_dict[original_tax_id]
	#			context.accepted_species_names_by_synonym_id[original_tax_id]

	Spc_Descendants = context.UserFlags_dict['SpeciesDescendants']
	suspected_as_HighRanked=dict()
	list_high_ranked=list()

	names_taxid_f=open(context.species_names_taxid_file,mode="w",encoding='utf8',newline='')
	writer = csv.writer(names_taxid_f,delimiter=',')
	writer.writerow(['TaxID', 'Species_Name'])

	# Decide which DataBase to use for Name matching & Name resolution
	if context.UserFlags_dict['NR_DB_name'] == 'NCBI': #default
		db_name = '/groups/itay_mayrose/share/ploidb/Ploidb_CODE/shared_files/NCBI_Names/NCBI_Names_for_NR_July2017.db'
	elif context.UserFlags_dict['NR_DB_name'] == 'PLT': # the plant list
		db_name = ploidb_config['general']['TPL_ALL_DB']
	else:
		logger.debug("Wrong DataBase for Name resolution")
		sys.exit()
		return


	TaxID_Name_dict=dict()
	TaxID_input_list=[]
	names_for_nr=[]
	Names_input_list=[]
	Names_AfterNR_list=[]
	Request_TaxIds_list=[]
	INPUT_NameTax_dict=dict()

	# 1. Split input list into 2 lists: Names and TaxIds
	for item in input_list:
		if item.isdigit():
			TaxID_input_list.append(int(item))
		else:
			Names_input_list.append(item)


	# 2. NAME RESOLUTION section - TBD TBD !!!!!!!!!!!!!!!!!
	#----------------------------------------------------------
	#Handle the names according to user request, only relevant for species names (more than 1 str):
	if context.UserFlags_dict['NameResType'] == 'SpellingCorrection' or context.UserFlags_dict['NameResType'] == 'MatchedNames':
		for name in Names_input_list:
			if ' ' in name:
				names_for_nr.append(name)
		Names_input_list = list(set(Names_input_list) - set(names_for_nr))
		Names_AfterNR_list = db_name_matching(db_name,names_for_nr,context)
		print(Names_AfterNR_list)
	#if context.UserFlags_dict['NameResType'] == 'MatchedNames':
		# Also prepare Syn/Accepted dictionaries:
		#Check_for_Synonyms()
		Names_input_list.extend(Names_AfterNR_list)


	logger.debug(Names_input_list)

	logger.debug("CheckPoint")
	logger.debug("CheckPoint")
	logger.debug("CheckPoint")
	# 3. Get TaxId for Names - input names Dictionary:
	NamesSTRING_for_query = turn_list_to_DB_list(Names_input_list)
	query_for_Names = "select * from ncbi_names_db where lower(Name) IN (%s)" %NamesSTRING_for_query
	db_rows = query_ncbi_db(query_for_Names,db_name)
	TaxId_list,Names_list,TaxId_rank_dict,Name_type_dict,TaxID_Name_dict,Name_TaxId_dict,TaxId_Parent_dict = extract_names_from_rows(context,db_rows,db_name)
	logger.debug("Names Vs TaxId(may be duplicates of TaxId, different names)")
	logger.debug(TaxID_Name_dict)
	logger.debug(Name_TaxId_dict)
	#Look for TaxId's with more than one name:
	list_TaxIds_manyVals=[]
	duple_taxidNames = {}
	for key, value in Name_TaxId_dict.items():
		duple_taxidNames.setdefault(value, set()).add(key)
	list_TaxIds_manyVals = [values for key, values in duple_taxidNames.items() if len(values) > 1]
	logger.debug(list_TaxIds_manyVals)
	#create the input dictionary:
	for name in Name_TaxId_dict.keys():
		logger.debug("----> name is %s" %name)
		TaxID = Name_TaxId_dict[name]
		logger.debug("Input name TaxId: %s, %s" %(name,TaxID))
		if name in Names_input_list:
			sequences_path = getSpeciesSequenceFileLocation(TaxID)  #check path for TaxId sequences in Genbank
			if sequences_path is None:
				# This can happen in cases where the taxon itself doesn't have sequences but its subspecies do
				logger.warning("Our Genbank has no sequences files for Tax ID " + str(TaxID) + "(Either no data or High Rank Taxon).")
				if ' ' in name: # Assuming High rank is with no space in name
					Names_input_list.remove(name)
			else:
				#check for more names with the same TaxId so we'll remove all and use the scientific name:
				if name in list_TaxIds_manyVals:
					logger.debug("Check name %s since its TaxId apears with more names in the list" %name)
					logger.debug("Type is %s" %Name_type_dict[name])
					if Name_type_dict[name] == 'scientific name':
						# Add check for scientific name so we'll keep it in case of few names input !!!!!!!!!!!!!!!
						context.species_names_by_tax_id[TaxID] = name
						if Spc_Descendants == 'No':
							Names_input_list.remove(name)
				else:
					context.species_names_by_tax_id[TaxID] = name
					if Spc_Descendants == 'No':
						Names_input_list.remove(name)

	logger.debug("Names_input_list")
	logger.debug(Names_input_list)
	logger.debug("context.species_names_by_tax_id")
	logger.debug(context.species_names_by_tax_id)






	###for name in Names_input_list:
	###	logger.debug("----> name is %s" %name)
	###	if name in Name_TaxId_dict.keys():  #Check if this name exists, if not remove
	###		TaxID = Name_TaxId_dict[name]
	###		INPUT_NameTax_dict[name] = TaxID
	###		logger.debug("Input name TaxId: %s, %s" %(name,TaxID))
	###		sequences_path = getSpeciesSequenceFileLocation(TaxID)  #check path for TaxId sequences in Genbank
	###		if sequences_path is None:
	###			# This can happen in cases where the taxon itself doesn't have sequences but its subspecies do
	###			logger.warning("Our Genbank has no sequences files for Tax ID " + str(TaxID) + "(Either no data or High Rank Taxon).")
	###		else:
	###			#check for more names with the same TaxId so we'll remove all and use the scientific name:
	###			Names_input_list.remove(name)
	###			context.species_names_by_tax_id[TaxID] = name
	###	else:
	###		logger.debug("Name %s does not exists in chosen DataBase and so it will be removed from the input list" %name)
	###		Names_input_list.remove(name)
###
	###logger.debug(Names_input_list)
	###logger.debug(context.species_names_by_tax_id)



	# 3. Convert Names to TaxIDs:
	NamesSTRING_for_query = turn_list_to_DB_list(Names_input_list)
	#query_for_Names = "select * from ncbi_names_db where Name IN (%s) and Type is 'scientific name'" %NamesSTRING_for_query
	query_for_Names = "select * from ncbi_names_db where lower(Name) IN (%s)" %NamesSTRING_for_query
	db_rows = query_ncbi_db(query_for_Names,db_name)
	TaxId_list,Names_list,TaxId_rank_dict,Name_type_dict,TaxID_Name_dict,Name_TaxId_dict,TaxId_Parent_dict = extract_names_from_rows(context,db_rows,db_name)
	#TaxId_list,Names_list,TaxId_rank_dict,TaxID_Name_dict,TaxId_Parent_dict = extract_names_from_rows_debug(context,db_rows,db_name)
	if TaxID_input_list:
		Request_TaxIds_list = set(TaxID_input_list + TaxId_list)
	else:
		Request_TaxIds_list=list(set(TaxId_list))
	logger.debug("User requested TaxIds list:")
	logger.debug(Request_TaxIds_list)


	#combine Name HR with TaxId HR
	#Working with Higher ranks
	#by distinguishing between parents with spc kids and parents with higher ranked kids:
	TaxID_input_list=[]
	#logger.debug("Initial list -> Spc List - TaxID_input_list:")
	#logger.debug(TaxID_input_list)


	if Spc_Descendants == 'Yes':
		TaxId_list_query=turn_list_to_DB_list(Request_TaxIds_list)
		#print(TaxId_list_query)
		#query = "select * from ncbi_names_db where ParentTaxID IN (%s) and Type is 'scientific name'" %TaxId_list_query
		query = "select * from ncbi_names_db where TaxID IN (%s) and Type is 'scientific name'" %TaxId_list_query
		db_rows = query_ncbi_db(query,db_name)
		TaxID_input_list,Names_list,TaxId_rank_dict,TaxID_Name_dict,Name_TaxId_dict,TaxId_Parent_dict = extract_names_from_rows_debug(context,db_rows,db_name)
		for item in TaxID_Name_dict.keys():
			logger.debug(item)
			logger.debug(TaxID_Name_dict[item])
		#sys.exit()

		# Split list into 2: species and higher rank.
		for TaxId in TaxID_Name_dict.keys():
			if ' ' not in TaxID_Name_dict[TaxId]:# and TaxId not in TaxID_Name_dict.keys():
				logger.debug("%s,-%s\n" % (str(TaxId),TaxID_Name_dict[TaxId]))
				suspected_as_HighRanked[TaxId] = (TaxID_Name_dict[TaxId])
				list_high_ranked.append(TaxId)
		#get all kids for high ranked TaxIds:
		logger.debug("Get all the descendants of all high ranked TaxIds:")
		while list_high_ranked:
			temp_HR_list = list_high_ranked
			Spc_rank_TaxId,High_rank_TaxId = Split_low_high_Ranks_debug(context,temp_HR_list,db_name)
			logger.debug("Spc List - Spc_rank_TaxId:")
			logger.debug(Spc_rank_TaxId)
			logger.debug("High rank taxa, need to check for kids:")
			logger.debug(High_rank_TaxId)
			#Update theTaxId species list:
			for spc_to_add in Spc_rank_TaxId:
				if spc_to_add not in TaxID_input_list:
					TaxID_input_list.append(spc_to_add)
			#Continue with parents taxIds:
			list_high_ranked = list(High_rank_TaxId)
	else:
		while Request_TaxIds_list:
			temp_HR_list = Request_TaxIds_list
			Spc_rank_TaxId,High_rank_TaxId = Split_low_high_Ranks_debug(context,temp_HR_list,db_name)
			logger.debug("Spc List - Spc_rank_TaxId:")
			logger.debug(Spc_rank_TaxId)
			logger.debug("Spc List - TaxID_input_list:")
			logger.debug(TaxID_input_list)
			logger.debug("High rank taxa, need to check for kids:")
			logger.debug(High_rank_TaxId)
			#Update theTaxId species list:
			for spc_to_add in Spc_rank_TaxId:
				if spc_to_add not in TaxID_input_list:
					TaxID_input_list.append(spc_to_add)
			#Continue with parents taxIds:
			Request_TaxIds_list = list(High_rank_TaxId)


	#When done
	logger.debug("TaxId list - all")
	logger.debug(TaxID_input_list)
	TaxID_input_list_str = turn_list_to_DB_list((TaxID_input_list))
	query = "select * from ncbi_names_db where TaxID IN (%s) and Type is 'scientific name'" %TaxID_input_list_str
	db_rows = query_ncbi_db(query,db_name)
	TaxId_list,Names_list,TaxId_rank_dict,Name_type_dict,TaxID_Name_dict,Name_TaxId_dict,TaxId_Parent_dict = extract_names_from_rows(context,db_rows,db_name)
	with open(context.working_dir+'/input_list_ncbi_data.csv','w') as input_data_ncbi:
		for line in db_rows:
			input_data_ncbi.write('"'+str(line[0])+'",'+'"'+str(line[1])+'",'+'"'+str(line[2])+'",'+'"'+str(line[3])+'",'+'"'+str(line[5])+'"\n')
			context.Pipe_input_Spc_Names_vs_TaxId[str(line[3])] = str(line[1])





	#sys.exit()

	###if Spc_Descendants != 'Yes':
	###	#remove high ranked names from list
	###	TaxID_input_list_str = turn_list_to_DB_list((TaxID_input_list))
	###	#Set query according to user descendent level:
	###	#if context.UserFlags_dict['DescendantLevel'] != 'specieLevel':
	###	query = "select * from ncbi_names_db where TaxID IN (%s) and Type is 'scientific name' and rank IN ('species','subspecies','varietas')" %TaxID_input_list_str
	###	#else:
	###	#	query = "select * from ncbi_names_db where TaxID IN (%s) and Type is 'scientific name' and rank IN ('species')" %TaxID_input_list_str
	###	db_rows = query_ncbi_db(query,db_name)
	###	TaxId_list,Names_list,TaxId_rank_dict,Name_type_dict,TaxID_Name_dict,Name_TaxId_dict,TaxId_Parent_dict = extract_names_from_rows(context,db_rows,db_name)
###
	###	logger.debug("No High rank in list")
	###	logger.debug(TaxId_list)
	###	#sys.exit()
###
	###	#get sientific names for TaxIds:
	###	#keys_list=[]
	###	#for key in TaxId_rank_dict.keys():
	###	#	keys_list.append(key)
	###	#keys_list_for=turn_list_to_DB_list(keys_list)
	###	tax_scient_str = turn_list_to_DB_list(TaxId_list)
	###	query = "select * from ncbi_names_db where TaxID IN (%s) and Type is 'scientific name'" %tax_scient_str
	###	db_rows = query_ncbi_db(query,db_name)
	###	TaxId_list,Names_list,TaxId_rank_dict,Name_type_dict,TaxID_Name_dict,Name_TaxId_dict,TaxId_Parent_dict = extract_names_from_rows(context,db_rows,db_name)

	#For debug Only:
	with open(context.working_dir+'/outNR.csv','w') as f_out:
		for key in TaxId_rank_dict.keys():
			#f_out.write('%s,%s,%s\n' %(str(key),str(TaxID_Name_dict[key]).encode("latin-1"),str(TaxId_rank_dict[key])))
			f_out.write('%s,%s,%s\n' %(str(key),str(TaxID_Name_dict[key]),str(TaxId_rank_dict[key])))
			context.species_names_by_tax_id[str(key)] = str(TaxID_Name_dict[key])
			writer.writerow([str(key), str(TaxID_Name_dict[key])])
	f_out.close()

	#Set the lists for OTT pipeline:
	context.species_list_ids = TaxId_list[:]
	context.species_list_names = Names_list[:]
	#context.parent_species_tax_id_by_subsp_tax_id
	#context.parent_species_names_by_subsp_tax_id
	logger.debug("NEW NEW NEW:")
	logger.debug(context.species_list_ids)
	logger.debug(context.species_list_names)
	#Once we have the TaxID list we need to extract all subsp and their parents for merging:
	logger.debug("SYNONYM:")
	logger.debug(context.syn_acc_TaxID_dict)
	logger.debug(context.syn_acc_dict)


	#Filter Results:
	temp_dict = copy.deepcopy(context.Pipe_input_Spc_Names_vs_TaxId)
	for name in temp_dict.keys():
		if context.UserFlags_dict['Filter_SubSp'] == 'on':
			logger.debug("Filter_SubSp is on")
			if ' subsp ' in name or ' var ' in name:
				#context.Pipe_input_Spc_Names_vs_TaxId.pop(name, None)
				del context.Pipe_input_Spc_Names_vs_TaxId[name]
		if context.UserFlags_dict['Filter_Hybrids'] == 'on':
			logger.debug(name)
			if ' x ' in name:
				del context.Pipe_input_Spc_Names_vs_TaxId[name]
	#since we have names with both 'x' and 'cf' we need to recopy the dictionary to exclude the names already removed.
	temp_dict = copy.deepcopy(context.Pipe_input_Spc_Names_vs_TaxId)
	for name in temp_dict.keys():
		if context.UserFlags_dict['Filter_Unresolved'] == 'on':
			if ' cf. ' in name or ' sp. ' in name or ' aff. ' in name or ' f. ' in name:
				del context.Pipe_input_Spc_Names_vs_TaxId[name]

	logger.debug("After Filter: Pipe_input_Spc_Names_vs_TaxId")
	logger.debug(context.Pipe_input_Spc_Names_vs_TaxId)
	sys.exit()

	return


#This function will create an output file for
def init_syn_accepted_Dict(context):

	syn_acc_dict=dict()
	syn_acc_TaxID_dict=dict()
	f_out = open('/groups/itay_mayrose/share/ploidb/ploidb_DBs/Syn_Acc_TaxId.csv','w')
	f_out.write('Synonym,Accepted\n')
	counter=0
	ncbi = NCBITaxa()
	#with open (f, 'r', newline='', encoding='latin-1') as infile
	df=pd.read_csv('/groups/itay_mayrose/share/ploidb/ploidb_DBs/Synonym_Table_PlantList.csv',encoding='latin-1')
	syn_names = df['Specie']
	acc_names = df['Accepted_Specie']
	len_data = len(syn_names)
	for ind in range(0,len_data):
		counter+=1
		syn_acc_dict[syn_names[ind]] = acc_names[ind]
		#Get TaxId for syn:
		syn_item = syn_names[ind].rstrip()
		nametranslator = ncbi.get_name_translator([syn_item])
		if nametranslator:
			syn_id = ncbi.get_name_translator([syn_item])[syn_item][0]
		else:
			#print("No TaxId for synonym name %s" %syn_item)
			continue
		#Get TaxId for syn:
		acc_item = acc_names[ind].rstrip()
		nametranslator = ncbi.get_name_translator([acc_item])
		if nametranslator:
			acc_id = ncbi.get_name_translator([acc_item])[acc_item][0]
		else:
			#print("No TaxId for accepted name %s" %acc_item)
			continue
		syn_acc_TaxID_dict[syn_id] = acc_id
		f_out.write(str(syn_id) +',' + str(acc_id)+'\n')
		print(counter)


def get_Accepted_name(original_name):

	db_Syn_Acc = ploidb_config['name resolve']['syn_acc_tpl']

	#db_NCBI_nodes = ploidb_config['general']['NCBI_Nodes_DB']
	conn_SynAcc_db = sqlite3.connect(db_Syn_Acc)
	curs_SynAcc = conn_SynAcc_db.cursor()

	query_getAccepted = "SELECT Accepted from Syn_Acc_table where Synonym like '%s'" %original_name
	curs_SynAcc.execute(query_getAccepted)
	AcceptedName_rows_db = curs_SynAcc.fetchall()

	if AcceptedName_rows_db:
		for item in AcceptedName_rows_db:
			Accepted_name = str(item[0])
		logger.debug("Accepted name for original '%s' is '%s'" %(original_name,Accepted_name))
	else:
		Accepted_name = 'None'

	return Accepted_name


#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
# NOT USED -> Rerun : Add remove species in case user edited taxa list
def Add_Remove_Species_from_clusters(context):

	#Get sequences for New species:
	if len(context.rerun_add_species_list) > 0: #in case we need to add species:
		create_dir_if_not_exists(context.rerun_new_taxIds_dir)
		create_dir_if_not_exists(context.rerun_new_taxIds_dir + '/Without_ITS/')
		create_dir_if_not_exists(context.rerun_new_taxIds_dir + '/ITS_only/')

		rerun_NewSpecies_Ids, rerun_NewSpecies_names, parents_list, parents_byIds = \
		get_species_ids_and_names_from_taxa_list(context,context.rerun_add_species_list,False,False)
		logger.info("New species List:")
		logger.info(rerun_NewSpecies_names)
		logger.info("New species Ids:")
		logger.info(rerun_NewSpecies_Ids)

		p = PathHelper(context.working_dir, context.id)
		f_log_list = open(context.rerun_new_taxIds_dir + '/LogList_AddedTaxIds.txt', 'w')
		#f_added_TaxId_gis = open(context.summary_dir + '/Added_NewTaxIds.txt', 'w')
		f_added_TaxId_gis = open(context.rerun_new_taxIds_dir + '/Added_NewTaxIds.txt', 'w')
		for taxId in rerun_NewSpecies_Ids:
			withITS_flag = False
			withoutITS = False
			fasta_newTax_filename = context.rerun_new_taxIds_dir + '/New_TaxId_' + str(taxId) + '_seqs.fasta'
			fasta_newTax_ITS_f = context.rerun_new_taxIds_dir + '/New_TaxId_' + str(taxId) + '_ITS.fasta'
			fasta_newTax_WithoutITS_f = context.rerun_new_taxIds_dir + '/New_TaxId_' + str(taxId) + '_NoITS.fasta'
			Dict_Add_gi_to_Cluster={} #keys are clusters, values are gi's to add
			with open(context.rerun_new_taxIds_dir+ "/LOG_blastall_" + str(taxId) + ".txt",'w') as log_file:
				log_file.write("ClusterNum,BitScore_cutOff,BestBitScore_NewSpecie,NewSpecie_gi\n")
				list_taxIds = [taxId]
				get_taxa_sequences(context,list_taxIds, fasta_filename = fasta_newTax_filename, gb_filename=context.gb_new_seq_filename, p=p)
				f_log_list.write(fasta_newTax_filename)
				f_log_list.write('\n')
				#Remove ITS seqs - handle later:
				remove_its_seqs(context,fasta_newTax_filename,fasta_newTax_ITS_f,fasta_newTax_WithoutITS_f)
				illegal_seqs=[]
				if os.path.exists(fasta_newTax_WithoutITS_f):
					withITS_flag = True
					#Run formatdb:
					if not os.path.exists(fasta_newTax_WithoutITS_f+'.nhr'):
						formatdb_cmd = "formatdb -i " + fasta_newTax_WithoutITS_f + " -pF -o T"
						os.system(formatdb_cmd)
				if os.path.exists(fasta_newTax_ITS_f):
					withoutITS = True
					if not os.path.exists(fasta_newTax_ITS_f+'.nhr'):
						formatdb_cmd = "formatdb -i " + fasta_newTax_ITS_f + " -pF -o T"
						os.system(formatdb_cmd)
				logger.debug("For TaxID %s: WithoutITS=%s, WithITS=%s !!!!!!!!!!!!!!!!" %(str(taxId),withoutITS,withITS_flag))

				#fasta_add_filename = context.rerun_new_taxIds_dir + '/Added_species_seqs.fasta'
				for cluster_id in context.rerun_chosen_clusters:
					path_of_selected_seq = context.working_dir + "/" + str(cluster_id) + "/seq-to-blast"
					with open(path_of_selected_seq,'r') as f_seq:
						logger.debug("Get data of selected seq of cluster %s" %str(cluster_id))
						first_line_seq = f_seq.readline()
						start_seqid='seqid|'
						end_seqid='|description|'
						seq_id = first_line_seq[first_line_seq.find(start_seqid)+len(start_seqid):first_line_seq.rfind(end_seqid)]
						start_gi = '>gi|'
						end_gi = '|taxonid|'
						gi_num = first_line_seq[first_line_seq.find(start_gi)+len(start_gi):first_line_seq.rfind(end_gi)]
						logger.debug("Selected seq: gi %s, seqid %s" %(gi_num,seq_id))

					#file_of_NewSpecie = context.rerun_new_taxIds_dir + '/New_species_' + str(taxId) + '_seqs.fasta'
					blast_results = context.rerun_new_taxIds_dir + "/seq_vs_NewSeq_" + str(taxId) + "_clst_" + str(cluster_id) + ".blastn"
					if os.path.exists(context.working_dir + "/" + str(cluster_id) + '/ITS_CLUSTER'):
						logger.debug("Cluster %s is an ITS cluster, use ITS seqs only" %str(cluster_id))
						if Add_ITS_seqs(context) and withITS_flag:
							blast_cmd = "blastall -p blastn -d " + fasta_newTax_ITS_f + " -i " + path_of_selected_seq + " -v 150000 -b 150000 -e 0.1 -m 8 > " + blast_results   # -m 9 will give a header to the output file
						else:
							continue
					else:
						if withoutITS:
							blast_cmd = "blastall -p blastn -d " + fasta_newTax_WithoutITS_f + " -i " + path_of_selected_seq + " -v 150000 -b 150000 -e 0.1 -m 8 > " + blast_results   # -m 9 will give a header to the output file
						else:
							continue
					logger.debug("Running blast cmd: %s" % blast_cmd)
					os.system(blast_cmd)
					blast_FilterOutput = context.rerun_new_taxIds_dir + "/seq_vs_NewSeq_" + str(taxId) + "_clst_" + str(cluster_id) + "filtered.blastn"
					log_file_filter = context.rerun_new_taxIds_dir + "/log_filter_"+ str(taxId) + "_clst" + str(cluster_id)+".txt"
					remove_shortLong_seqs(blast_results, fasta_newTax_WithoutITS_f, path_of_selected_seq, blast_FilterOutput,log_file_filter)
					cutOff_BitScore = check_cluster_BlastAll(cluster_id,seq_id,context.working_dir)
					if os.stat(blast_FilterOutput).st_size == 0 or cutOff_BitScore == 999:
						largeSeq_blastResults = [0,0]
					else:
						largeSeq_blastResults = checkLargeSeqBlast(blast_FilterOutput)
					#Add seq to aligned file:
					log_file.write(str(cluster_id) + ',' + str(cutOff_BitScore)+ ',' + str(largeSeq_blastResults[1]) + ',' + str(largeSeq_blastResults[0])+'\n')
					#Check if BitScore is >= to CutOff value:
					if (float(largeSeq_blastResults[1]) >= cutOff_BitScore):
						f_added_TaxId_gis.write("\nTaxID-%s: " % str(taxId))
						if largeSeq_blastResults[0] in Dict_Add_gi_to_Cluster.values():
							illegal_seqs.append(largeSeq_blastResults[0])
						else:
							Dict_Add_gi_to_Cluster[cluster_id] = largeSeq_blastResults[0]
							f_added_TaxId_gis.write("%s - %s," %(str(cluster_id),str(largeSeq_blastResults[0])))
							gi_num = largeSeq_blastResults[0].split('|')[1]
							#Create seq file for adding:
							add_seq_file = context.rerun_new_taxIds_dir + "/add_taxId_" + str(taxId) + "_gi_" + str(gi_num) + "_toClust_"+ str(cluster_id) + '.fasta'
							add_seq_withDesc_file = context.rerun_new_taxIds_dir + "/add_desc_taxId_" + str(taxId) + "_gi_" + str(gi_num) + "_toClust_"+ str(cluster_id) + '.fasta'
							with open(fasta_newTax_WithoutITS_f, "rU") as fasta_newTaxId_f:
								all_seqs = list(SeqIO.parse(fasta_newTaxId_f, "fasta"))
								for seq in all_seqs:
									gi_num_identify = 'gi|'+gi_num
									if gi_num_identify in seq.description:
										with open(add_seq_file,'w') as f_seqToAdd:
											SeqIO.write(seq, f_seqToAdd, "fasta")
											#organism_name = (getPropertyFromFastaSeqHeader(seq.description, "organism")).replace(' ','_')
											organism_name = return_codedName(getPropertyFromFastaSeqHeader(seq.description, "organism"))
											f_seqToAdd.close()
										shutil.copyfile(add_seq_file, add_seq_withDesc_file)
										replaceToOrganismName(add_seq_file,organism_name)
						fasta_newTaxId_f.close()

				#Check each cluster was selected once:
				logger.debug(Dict_Add_gi_to_Cluster)
				if not Dict_Add_gi_to_Cluster:
					logger.debug("No matched clusters for Added Species")
					with open(context.summary_file,'a') as f_sum:
						f_sum.write("None of the New taxa matched the data in current clusters\n")
						#return
				else:
					logger.debug("New taxa gi's added to clusters - Dict_Add_gi_to_Cluster:")
					logger.debug(Dict_Add_gi_to_Cluster)

				if illegal_seqs:
					logger.debug("New taxID matched more than one cluster, need to verify different gi's !!")
					logger.debug("Illegal seqs:")
					logger.debug(illegal_seqs)
					return
				else:
					for cluster_id in context.rerun_chosen_clusters:
						if cluster_id in Dict_Add_gi_to_Cluster.keys():
							gi_num = Dict_Add_gi_to_Cluster[cluster_id].split('|')[1]
							add_seq_file = context.rerun_new_taxIds_dir + "/add_taxId_" + str(taxId) + "_gi_" + str(gi_num) + "_toClust_"+ str(cluster_id) + '.fasta'
							add_seq_withDesc_file = context.rerun_new_taxIds_dir + "/add_desc_taxId_" + str(taxId) + "_gi_" + str(gi_num) + "_toClust_"+ str(cluster_id) + '.fasta'
							#Add seq to existing all seqs file:
							seqs_noMulti_cluster_f = context.working_dir + '/' +str(cluster_id) + "/seqs-no-multiple-accessions.fasta"
							shutil.copyfile(seqs_noMulti_cluster_f,seqs_noMulti_cluster_f +'_original')
							logger.debug("Merge seqs files: %s and %s" %(seqs_noMulti_cluster_f,seqs_noMulti_cluster_f))
							merge_seq_files([seqs_noMulti_cluster_f, add_seq_withDesc_file],seqs_noMulti_cluster_f)
							replaceToLowerCase(seqs_noMulti_cluster_f)

	#Remove species from clusters:
	#for cluster_id in context.rerun_chosen_clusters:
	#	for specie_toRemove in context.rerun_remove_species_list:
	#		#fasta_file_toEdit = '/bioseq/data/results/oneTwoTree/1486037489/concat/2/seqs-with-out-group-concat.fasta'
	#		fasta_file_toEdit = context.working_dir + '/concat/' + str(cluster_id) + '/seqs-with-out-group-concat.fasta'
	#		remove_seq_from_file(specie_toRemove, fasta_file_toEdit)


	return

#-----------------------------------------------------------------------------------------------------------
def Split_ITSseqs_NewTaxIds(context):
#This function will take all new taxIds and collect its seqs. Once we have that we'll split them into ITS and NonITS
# fasta files to handle their addition to the existing clusters

	#Get sequences for New species:
	if len(context.rerun_add_species_list) > 0 and context.rerun_add_species_list[0] is not '': #in case we need to add species:
		create_dir_if_not_exists(context.rerun_new_taxIds_dir)
		NonITS_dir = context.rerun_new_taxIds_dir + '/Without_ITS/'
		ITS_dir = context.rerun_new_taxIds_dir + '/ITS_only/'
		create_dir_if_not_exists(ITS_dir)
		create_dir_if_not_exists(NonITS_dir)

		rerun_NewSpecies_Ids, rerun_NewSpecies_names, parents_list, parents_byIds = \
		get_species_ids_and_names_from_taxa_list(context,context.rerun_add_species_list,False,False)
		logger.info("New species List:")
		logger.info(rerun_NewSpecies_names)
		logger.info("New species Ids:")
		logger.info(rerun_NewSpecies_Ids)

		p = PathHelper(context.working_dir, context.id)
		f_log_list = open(context.rerun_new_taxIds_dir + '/AddedTaxIds_filePath.txt', 'w')
		illegal_seqs=[]
		for taxId in rerun_NewSpecies_Ids:
			Dict_Add_gi_to_Cluster={}   			#keys are clusters, values are gi's to add
			#seqs files names:
			fasta_newTax_filename = context.rerun_new_taxIds_dir + '/New_TaxId_' + str(taxId) + '_seqs.fasta'
			fasta_newTax_ITS_f = ITS_dir + '/New_TaxId_' + str(taxId) + '_ITS.fasta'
			fasta_newTax_WithoutITS_f = NonITS_dir + '/New_TaxId_' + str(taxId) + '_NoITS.fasta'
			list_taxIds = [taxId]
			get_taxa_sequences(context,list_taxIds, fasta_filename = fasta_newTax_filename, gb_filename=context.gb_new_seq_filename, p=p)
			f_log_list.write(fasta_newTax_filename)
			f_log_list.write('\n')
			#Remove ITS seqs and run formatdb:
			remove_its_seqs(context,fasta_newTax_filename,fasta_newTax_ITS_f,fasta_newTax_WithoutITS_f)
			#The file to which we'll write the results of each blast:
			log_file = open(context.rerun_new_taxIds_dir+ "/LOG_blastall_" + str(taxId) + ".txt",'w')
			log_file.write("ClusterNum,BitScore_cutOff,BestBitScore_NewSpecie,NewSpecie_gi\n")

			#Go through all selected clusters to check for a match:
			for cluster_id in context.rerun_chosen_clusters:
				path_of_selected_seq = context.working_dir + "/" + str(cluster_id) + "/seq-to-blast"
				with open(path_of_selected_seq,'r') as f_seq:
					logger.debug("Get data of selected seq of cluster %s" %str(cluster_id))
					first_line_seq = f_seq.readline()
					start_seqid='seqid|'
					end_seqid='|description|'
					seq_id = first_line_seq[first_line_seq.find(start_seqid)+len(start_seqid):first_line_seq.rfind(end_seqid)]
					start_gi = '>gi|'
					end_gi = '|taxonid|'
					gi_num = first_line_seq[first_line_seq.find(start_gi)+len(start_gi):first_line_seq.rfind(end_gi)]
					logger.debug("Selected seq: gi %s, seqid %s" %(gi_num,seq_id))


				if os.path.exists(context.working_dir + "/" + str(cluster_id) + '/ITS_CLUSTER'):
					logger.debug("Cluster %s is an ITS cluster, use ITS seqs only" %str(cluster_id))
					#Check if there's seqs file for this taxid:
					if os.path.exists(fasta_newTax_ITS_f):
						#AddFragments
						blast_results = ITS_dir + "/blastall_TaxId_" + str(taxId) + "_clstITS_" + str(cluster_id) + ".blastn"
						blast_cmd = "blastall -p blastn -d " + fasta_newTax_ITS_f + " -i " + path_of_selected_seq + " -v 150000 -b 150000 -e 0.1 -m 8 > " + blast_results   # -m 9 will give a header to the output file
						os.system(blast_cmd)
						#In case of more than 1 seq, need to handle multiple accessions:
						temp_seqs_list = list(SeqIO.parse(fasta_newTax_ITS_f, "fasta"))
						if len(temp_seqs_list) > 1:
							#run blast all-vs-all:
							blast_allAll_results_filename = fasta_newTax_ITS_f + '-allvsall'
							blast_command = "blastall -p blastn -d %s -i %s -v 100000 -b 100000 -e 1e-5 -m 8" % (
								fasta_newTax_ITS_f, fasta_newTax_ITS_f)
							exec_external_command_redirect_output(blast_command, blast_allAll_results_filename)
							#remove multiple accessions:
							rerun_get_representative_seq(taxId,fasta_newTax_ITS_f,fasta_newTax_ITS_f + '_rep',blast_allAll_results_filename)
						else:
							shutil.copyfile(fasta_newTax_ITS_f,fasta_newTax_ITS_f + '_rep')

						seqs_noMulti_cluster_f = context.working_dir + '/' +str(cluster_id) + "/seqs-no-multiple-accessions.fasta"
						shutil.copyfile(seqs_noMulti_cluster_f,seqs_noMulti_cluster_f +'_original')
						addFrag_cmd = 'mafft --addfragments ' + fasta_newTax_ITS_f + '_rep' + ' --multipair ' + seqs_noMulti_cluster_f + ' > ' + seqs_noMulti_cluster_f + '_added'
						os.system(addFrag_cmd)
						logger.debug("Taxon ID %s was added to cluster %s" %(taxId,cluster_id))
						shutil.copyfile(seqs_noMulti_cluster_f + '_added',seqs_noMulti_cluster_f)
						break
					else:
						logger.debug("No ITS seqs for TaxId %s........continue to next taxid" %str(taxId))
						break
				else:
					logger.debug("Cluster %s is a Non ITS cluster, use NonITS seqs only" %str(cluster_id))
					if os.path.exists(fasta_newTax_WithoutITS_f):
						blast_results = NonITS_dir + "/blastall_TaxId_" + str(taxId) + "_clstNonITS_" + str(cluster_id) + ".blastn"
						blast_cmd = "blastall -p blastn -d " + fasta_newTax_WithoutITS_f + " -i " + path_of_selected_seq + " -v 150000 -b 150000 -e 0.1 -m 8 > " + blast_results   # -m 9 will give a header to the output file
						blast_FilterOutput = NonITS_dir + "/blastall_TaxId_" + str(taxId) + "_clstNonITS_" + str(cluster_id) + "filtered.blastn"
						log_file_filter = NonITS_dir + "/log_filter_"+ str(taxId) + "_clstNonITS_" + str(cluster_id)+".txt"
						remove_shortLong_seqs(blast_results, fasta_newTax_WithoutITS_f, path_of_selected_seq, blast_FilterOutput,log_file_filter)
						fasta_seq_file = fasta_newTax_WithoutITS_f
					else:
						logger.debug("No NON-ITS seqs for TaxId %s....continue to next taxid" %str(taxId))
						break
				logger.debug("Running blast cmd: %s" % blast_cmd)

				os.system(blast_cmd)
				#blast_FilterOutput = context.rerun_new_taxIds_dir + "/seq_vs_NewSeq_" + str(taxId) + "_clst_" + str(cluster_id) + "filtered.blastn"
				#log_file_filter = context.rerun_new_taxIds_dir + "/log_filter_"+ str(taxId) + "_clst" + str(cluster_id)+".txt"
				#remove_shortLong_seqs(blast_results, fasta_newTax_WithoutITS_f, path_of_selected_seq, blast_FilterOutput,log_file_filter)
				cutOff_BitScore = check_cluster_BlastAll(cluster_id,seq_id,context.working_dir)
				#if os.stat(blast_FilterOutput).st_size == 0 or cutOff_BitScore == 999:
				if (not os.path.exists(blast_FilterOutput)) or (cutOff_BitScore == 999):
					largeSeq_blastResults = [0,0]
				else:
					largeSeq_blastResults = checkLargeSeqBlast(blast_FilterOutput)
				#Add seq to aligned file:
				log_file.write(str(cluster_id) + ',' + str(cutOff_BitScore)+ ',' + str(largeSeq_blastResults[1]) + ',' + str(largeSeq_blastResults[0])+'\n')
				#Check if BitScore is >= to CutOff value:
				if (float(largeSeq_blastResults[1]) >= cutOff_BitScore):
					f_added_TaxId_gis.write("\nTaxID-%s: " % str(taxId))
					if largeSeq_blastResults[0] in Dict_Add_gi_to_Cluster.values():
						illegal_seqs.append(largeSeq_blastResults[0])
					else:
						Dict_Add_gi_to_Cluster[cluster_id] = largeSeq_blastResults[0]
						f_added_TaxId_gis.write("%s - %s," %(str(cluster_id),str(largeSeq_blastResults[0])))
						gi_num = largeSeq_blastResults[0].split('|')[1]
						#Create seq file for adding:
						add_seq_file = context.rerun_new_taxIds_dir + "/add_taxId_" + str(taxId) + "_gi_" + str(gi_num) + "_toClust_"+ str(cluster_id) + '.fasta'
						add_seq_withDesc_file = context.rerun_new_taxIds_dir + "/add_desc_taxId_" + str(taxId) + "_gi_" + str(gi_num) + "_toClust_"+ str(cluster_id) + '.fasta'
						with open(fasta_seq_file, "rU") as fasta_newTaxId_f:
							all_seqs = list(SeqIO.parse(fasta_newTaxId_f, "fasta"))
							for seq in all_seqs:
								gi_num_identify = 'gi|'+gi_num
								if gi_num_identify in seq.description:
									with open(add_seq_file,'w') as f_seqToAdd:
										SeqIO.write(seq, f_seqToAdd, "fasta")
										#organism_name = (getPropertyFromFastaSeqHeader(seq.description, "organism")).replace(' ','_')
										organism_name = return_codedName(getPropertyFromFastaSeqHeader(seq.description, "organism"))
										f_seqToAdd.close()
									shutil.copyfile(add_seq_file, add_seq_withDesc_file)
									replaceToOrganismName(add_seq_file,organism_name)
					fasta_newTaxId_f.close()

		#Check each cluster was selected once:
		logger.debug(Dict_Add_gi_to_Cluster)
		if not Dict_Add_gi_to_Cluster:
			logger.debug("No matched clusters for Added Species")
			with open(context.summary_file,'a') as f_sum:
				f_sum.write("None of the New taxa matched the data in current clusters\n")
				#return
		else:
			logger.debug("New taxa gi's added to clusters - Dict_Add_gi_to_Cluster:")
			logger.debug(Dict_Add_gi_to_Cluster)

		if illegal_seqs:
			logger.debug("New taxID matched more than one cluster, need to verify different gi's !!")
			logger.debug("Illegal seqs:")
			logger.debug(illegal_seqs)
			return
		else:
			for cluster_id in context.rerun_chosen_clusters:
				if cluster_id in Dict_Add_gi_to_Cluster.keys():
					gi_num = Dict_Add_gi_to_Cluster[cluster_id].split('|')[1]
					add_seq_file = context.rerun_new_taxIds_dir + "/add_taxId_" + str(taxId) + "_gi_" + str(gi_num) + "_toClust_"+ str(cluster_id) + '.fasta'
					add_seq_withDesc_file = context.rerun_new_taxIds_dir + "/add_desc_taxId_" + str(taxId) + "_gi_" + str(gi_num) + "_toClust_"+ str(cluster_id) + '.fasta'
					#Add seq to existing all seqs file:
					seqs_noMulti_cluster_f = context.working_dir + '/' +str(cluster_id) + "/seqs-no-multiple-accessions.fasta"
					shutil.copyfile(seqs_noMulti_cluster_f,seqs_noMulti_cluster_f +'_original')
					logger.debug("Merge seqs files: %s and %s" %(seqs_noMulti_cluster_f,seqs_noMulti_cluster_f))
					merge_seq_files([seqs_noMulti_cluster_f, add_seq_withDesc_file],seqs_noMulti_cluster_f)
					replaceToLowerCase(seqs_noMulti_cluster_f)


		logger.debug("New ITS files to check:")
		for path in context.rerun_new_its_files_list:
			logger.debug(path)
		logger.debug("New NON-ITS files to check:")
		for path in context.rerun_new_non_its_files_list:
			logger.debug(path)


	return


# Michal Jan 2017
def get_species_ids_and_names_from_taxa_list(context,taxa_list, should_ignore_subsp=False, should_merge_subsp=False):
	logger.info("Converting taxa list to list of species ids and species names (%d)" %len(taxa_list))
	species_ids = []
	species_names = []
	parent_species_tax_id_by_subsp_tax_id = dict()
	parent_species_names_by_subsp_tax_id = dict()

	ncbi = NCBITaxa()
	#Check if digit -> TaxID taxon list:
	if taxa_list[0].isdigit():
		logger.info("Input taxa list was identifies as TaxID list")
		switch_toTaxId = 0
	else:
		logger.info("Input taxa list was identifies as Names list")
		switch_toTaxId = 1
	first_flag=0

	#get_species_TaxIds_and_Names(context,taxa_list, should_ignore_subsp=False, should_merge_subsp=False)
	#TaxId_list=[]   # All kids of all input TaxIds/Taxon names
	tax_id_dict={}
	for taxon_item in taxa_list:
		if not taxon_item.isdigit():
			taxon_item = taxon_item.rstrip()
			nametranslator = ncbi.get_name_translator([taxon_item])
			if nametranslator:
				taxon_id = ncbi.get_name_translator([taxon_item])[taxon_item][0]
				logger.debug(taxon_id)
				logger.debug(ncbi.get_name_translator([taxon_item]))
			else:
				logger.debug("No TaxId for taxon %s" %taxon_item)
				continue
		else:
			taxon_id = taxon_item
		taxon_rank = ncbi.get_rank([int(taxon_id)])#  -> {143: 'species'}
		logger.info(taxon_id)
		logger.info(taxon_rank)
		cur_ids = get_allKids_forTaxId(str(taxon_id),context.UserFlags_dict['Include_SubSp'])
		if not cur_ids:
			logger.debug('No kids for TaxId %s: add taxId of parent %s' %(taxon_id,taxon_rank))
			species_ids.append(str(taxon_id))
			taxon_name = ncbi.get_taxid_translator([int(taxon_id)])[int(taxon_id)]
			tax_id_dict[str(taxon_id)] = taxon_name
			species_names.append(taxon_name)
		else:
			rank_parent = taxon_rank[int(taxon_id)]
			logger.debug('Found kids for TaxId %s (rank = %s): add taxIds of kids %s' %(taxon_id,rank_parent,cur_ids))
			#If parent is species than add it as well:
			if taxon_rank[int(taxon_id)] == 'species':
				species_ids.append(str(taxon_id))
				taxon_name = ncbi.get_taxid_translator([int(taxon_id)])[int(taxon_id)]
				tax_id_dict[str(taxon_id)] = taxon_name
				species_names.append(taxon_name)
			species_ids.extend(cur_ids)
			for taxon_id in cur_ids:
				taxon_name = ncbi.get_taxid_translator([int(taxon_id)])[int(taxon_id)]
				tax_id_dict[taxon_id] = taxon_name
				species_names.append(taxon_name)
		logger.debug("Current Species list:")
		logger.debug(species_names)
		logger.debug(species_ids)

	remove_tax_list=[]
	with open(context.species_names_taxid_file,mode="w",encoding='utf8',newline='') as names_taxid_f:
		writer = csv.writer(names_taxid_f,delimiter=',')
		writer.writerow(['TaxID', 'Species_Name'])
		#if first_flag == 0: writer.writerow(['TaxID', 'Species_Name'])
		for tax_id in species_ids:
			logger.debug("Checking NOW: %s" %tax_id)
			#Check if name is legal/ user asked to exclude (subsp. cf. and such):
			subsp_parent_name='None'
			retFlag_checkTaxa,subsp_parent_name = check_tax_id_name(context,tax_id_dict[tax_id],tax_id)
			if retFlag_checkTaxa is True:
				logger.debug("writeToFile: %s " %tax_id_dict[tax_id])
				writer.writerow([str(tax_id), str(tax_id_dict[tax_id])])
				context.species_names_by_tax_id[str(tax_id)] = tax_id_dict[tax_id]
				#In case there's subsp merge case:
				if subsp_parent_name != 'None':
					f_subsp_var_names=open(context.subsp_variants_merge_names,'a')
					parent_species_names_by_subsp_tax_id[tax_id]=subsp_parent_name
					parent_subVar_TaxId = (list(tax_id_dict.keys())[list(tax_id_dict.values()).index(subsp_parent_name)])
					parent_species_tax_id_by_subsp_tax_id[tax_id]=parent_subVar_TaxId
					f_subsp_var_names.write("%s,%s,%s,%s\n" % (tax_id,tax_id_dict[tax_id],parent_subVar_TaxId,subsp_parent_name))
			else:
				#remove this name from lists:
				remove_tax_list.append(tax_id)
				logger.debug("Remove name from list: %s, %s" % (tax_id_dict[tax_id],tax_id,))
		#remove illegal names according to user selection:
		for rem_tax in remove_tax_list:
			species_ids.remove(rem_tax)
			species_names.remove(tax_id_dict[rem_tax])

		#Add user outgroup as species in list:
		if context.UserFlags_dict['Outgroup_User'] == 'User':
			context.User_Outgroup = context.UserFlags_dict['Outgroup_User']
			#user_outgroup_taxid = Outgroup_User
			nametranslator = ncbi.get_name_translator([context.User_Outgroup])
			if nametranslator:
				user_outgroup_taxid = ncbi.get_name_translator([context.User_Outgroup])[context.User_Outgroup][0]
				logger.debug("User Outgroup Name: %s, TaxId: %s" %(context.User_Outgroup,user_outgroup_taxid))
				species_ids.append(user_outgroup_taxid)
				species_names.append(context.User_Outgroup)
				tax_id_dict[user_outgroup_taxid] = context.User_Outgroup
				logger.debug("writeToFile: %s " %tax_id_dict[user_outgroup_taxid])
				writer.writerow([str(user_outgroup_taxid), str(tax_id_dict[user_outgroup_taxid])])
				context.species_names_by_tax_id[str(user_outgroup_taxid)] = tax_id_dict[user_outgroup_taxid]
			else:
				with open(context.final_status,'a') as f_status:
					f_status.write("FAILED - User Outgroup TaxId wasn't found in NCBI, please check the Name\n")

	names_taxid_f.close()
	# clearing duplicates in the lists using set
	species_ids = list(set(species_ids))
	species_names = list(set(species_names))

	logger.info("Finished converting taxa list to list of species. Species list contain %d unique names" % len(
		species_names))

	return species_ids, species_names, parent_species_tax_id_by_subsp_tax_id, parent_species_names_by_subsp_tax_id



#Get names and taxIds of species according to the taxId given:
#1. extract the TaxId of the given name
#2. create list of descendants species using get_descendant_taxa ete3
#3. get species names according to TaxIds !!!
def get_species_from_taxon_name_IDver(context, taxon_name,first_flag,should_ignore_subsp=False, should_merge_subsp=False):

	names_taxid_f = open(context.species_names_taxid_file,mode="a",encoding='utf8',newline='')
	writer = csv.writer(names_taxid_f,delimiter=',')
	#Name,Authority,Original name,Original authority,Coded Name,Coded Authority,Score,Matched Name,Id
	#writer.writerow(['Name','Authority','Original name','Original authority','Coded Name','Coded Authority','Score','Matched Name','Id'])
	if first_flag==0: writer.writerow(['TaxID', 'Species_Name'])

	genera_flag=0
	species_names_by_tax_id = dict()
	parent_species_names_by_subsp_tax_id = dict()
	parent_species_tax_id_by_subsp_tax_id = dict()

	# TPL db name:
	db_name = ploidb_config['general']['NCBI_NAMES_TAXIDS_DB']  # Names of Sientific Names and TaxIds (from names.dmp)
	conn_tpl = sqlite3.connect(db_name)
	curs_tpl = conn_tpl.cursor()

	#1. extract TaxId of given name:
	# Check if name or digit (in case of TaxId input)
	if taxon_name.isdigit():
		main_tax_id = taxon_name
	else:
		query_taxId = "SELECT TaxID FROM names_vs_tax_ids WHERE Name LIKE '%s'" % taxon_name
		curs_tpl.execute(query_taxId)
		taxID_row = curs_tpl.fetchall()
		if taxID_row:
			for line in taxID_row:
				main_tax_id = str(line[0])
		else:
			logger.debug("Empty results for %s" %taxon_name)
			return list(), list(), dict(), dict()
	logger.debug("TaxId of taxon %s is %s" %(taxon_name,main_tax_id))


	# 2. create list of descendants species using get_descendant_taxa ete3
	ncbi = NCBITaxa()
	try:
		TaxId_list = set(ncbi.get_descendant_taxa(main_tax_id, collapse_subspecies=False,return_tree=False))
		rank_dict = ncbi.get_rank(main_tax_id)
		if len(TaxId_list) > 1:
			TaxId_list_nodes = set(list(ncbi.get_descendant_taxa(main_tax_id, collapse_subspecies=True,return_tree=False)))  #given in map
			diff_taxList  = TaxId_list_nodes - TaxId_list
			Final_TaxList = list(TaxId_list) + list(diff_taxList)
		else:
			Final_TaxList = list(TaxId_list)
		Final_TaxList.append(main_tax_id)
		Final_TaxList = set(Final_TaxList)
	except ValueError:
		Final_TaxList = [main_tax_id]
	logger.debug("TaxId list returned for %s: %s" %(taxon_name,Final_TaxList))

	query_line = "SELECT * FROM names_vs_tax_ids WHERE TaxID IN (%s)" % Final_TaxList
	query_line=query_line.replace('{','')
	query_line=query_line.replace('}','')
	query_line=query_line.replace('[','')
	query_line=query_line.replace(']','')
	logger.debug("query_line %s" %(query_line))
	query_name_syn = 'None'
	#query_line,query_name_syn, genera_flag=CreateTPLdb_query(taxon_name,context)

	curs_tpl.execute(query_line)
	rows_from_db = curs_tpl.fetchall()
	# Check if there's a syn query:
	No_syn=0
	if query_name_syn != 'None':
		No_syn=1
		curs_tpl.execute(query_name_syn)
		names_syn_fromdb = curs_tpl.fetchall()


	if('TRUE'):
		for line in rows_from_db:
			tax_id = str(line[0])			# TaxId of the original name (not of the accepted name)
			species_name = str(line[1])		# Accepted name - According to Name resolution (should match Coded name without the '_'
			#verify no cf. or sp. in species names:
			if 'cf.' in species_name or ' sp.' in species_name:
				logger.debug("Skip species %s (include either cf. or sp.)" % species_name)
			else:
				writer.writerow(line)
				context.species_names_by_tax_id[tax_id] = species_name
		#Check number of species:
		max_number_of_species = 20000
		number_of_species=len(context.species_names_by_tax_id)
		if number_of_species > max_number_of_species:
			logger.error("too many results for " + taxon_name + "(" + str(number_of_species) + "). Skipping")
			return list(), list()
		else:
			logger.debug("Number of species found in data base for Genus %s: %d" % (taxon_name,number_of_species))
			if number_of_species == 0:
				return list(), list(), dict(), dict()

		#Update Names and syn/name res file:
		if No_syn !=0:
			for line in names_syn_fromdb:
				name = str(line[0])
				accepted_name = str(line[1])
				#writer_syn.writerow([name,accepted_name])
				#context.names_vs_accepted_dict[name] = accepted_name


	names_taxid_f.close()
	#names_syn_f.close()

	#sys.exit(0)
	return list(context.species_names_by_tax_id.keys()), list(context.species_names_by_tax_id.values()), parent_species_tax_id_by_subsp_tax_id, parent_species_names_by_subsp_tax_id


#-----------------------------------
#--------------------------------------------------def-------------------------------------------------------
#This function will take the user's input list and return all names found in NCBI database
# until the descendant level defined by the user
def return_ncbi_species_list(input_list,context):
# Desc: this function will create the species list and all TaxId data according to user params:
# enables NCBI db alone, spelling correction and Acepted names matching

	# for Accepted names match need to prepare the correct dictionaries:
	# 1. for subsp merge:
	#			parent_tax_id = context.parent_species_tax_id_by_subsp_tax_id[original_tax_id]
	#			parent_name = context.parent_species_names_by_subsp_tax_id[original_tax_id]
	# 2. for Acc/Synonym merge:
	#			context.syn_acc_dict[original_tax_id]
	#			context.accepted_species_names_by_synonym_id[original_tax_id]

	#Get level of retrival:
	depth_level = 'species' #un the future can be either: 'species' (default) or 'subsp'

	names_taxid_f=open(context.species_names_taxid_file,mode="w",encoding='utf8',newline='')
	writer = csv.writer(names_taxid_f,delimiter=',')
	writer.writerow(['TaxID', 'Species_Name'])

	# Decide which DataBase to use for Name matching & Name resolution
	if context.UserFlags_dict['NR_DB_name'] == 'NCBI': #default
		db_name = '/groups/itay_mayrose/share/ploidb/Ploidb_CODE/shared_files/NCBI_Names/NCBI_Names_for_NR_July2017.db'
	elif context.UserFlags_dict['NR_DB_name'] == 'PLT': # the plant list
		db_name = ploidb_config['general']['TPL_ALL_DB']
	else:
		logger.debug("Wrong DataBase for Name resolution")
		sys.exit()
		return


	TaxID_Name_dict=dict()
	TaxID_input_list=[]
	names_for_nr=[]
	Names_input_list=[]
	Names_AfterNR_list=[]
	Request_TaxIds_list=[]

	# Check which inputs are taxIds and which Names
	for item in input_list:
		if item.isdigit():
			TaxID_input_list.append(int(item))
		else:
			Names_input_list.append(item)


	#NAME RESOLUTION section
	#----------------------------------------------------------
	#Handle the names according to user request, only relevant for species names (more than 1 str):
	if context.UserFlags_dict['NameResType'] == 'SpellingCorrection' or context.UserFlags_dict['NameResType'] == 'MatchedNames':
		for name in Names_input_list:
			if ' ' in name:
				names_for_nr.append(name)
		Names_input_list = list(set(Names_input_list) - set(names_for_nr))
		Names_AfterNR_list = db_name_matching(db_name,names_for_nr,context)
		print(Names_AfterNR_list)
	#if context.UserFlags_dict['NameResType'] == 'MatchedNames':
		# Also prepare Syn/Accepted dictionaries:
		#Check_for_Synonyms()
		Names_input_list.extend(Names_AfterNR_list)
	#sys.exit()

	# Convert Names to TaxIDs:
	NamesSTRING_for_query = turn_list_to_DB_list(Names_input_list)
	#query_for_Names = "select * from ncbi_names_db where Name IN (%s) and Type is 'scientific name'" %NamesSTRING_for_query
	query_for_Names = "select * from ncbi_names_db where lower(Name) IN (%s) and Type not in ('common name','authority')" %NamesSTRING_for_query
	db_rows = query_ncbi_db(query_for_Names,db_name)
	TaxId_list,Names_list,TaxId_rank_dict,Name_type_dict,TaxID_Name_dict,TaxId_Parent_dict = extract_names_from_rows(context,db_rows,db_name)
	if TaxID_input_list:
		Request_TaxIds_list = set(TaxID_input_list + TaxId_list)
	else:
		Request_TaxIds_list=list(TaxId_list)
	logger.debug("User requested TaxIds list:")
	logger.debug(Request_TaxIds_list)

	#combine Name HR with TaxId HR
	#Working with Higher ranks
	#by distinguishing between parents with spc kids and parents with higher ranked kids:
	TaxID_input_list=[]
	while Request_TaxIds_list:
		temp_HR_list = Request_TaxIds_list
		Spc_rank_TaxId,High_rank_TaxId = Split_low_high_Ranks(context,temp_HR_list,db_name)
		logger.debug("Spc List - Spc_rank_TaxId:")
		logger.debug(Spc_rank_TaxId)
		logger.debug("Spc List - TaxID_input_list:")
		logger.debug(TaxID_input_list)
		logger.debug("High rank taxa, need to check for kids:")
		logger.debug(High_rank_TaxId)
		#Update theTaxId species list:
		TaxID_input_list.extend(list(set(Spc_rank_TaxId)))
		#Continue with parents taxIds:
		Request_TaxIds_list = list(High_rank_TaxId)

	logger.debug("TaxId list - all")
	logger.debug(Request_TaxIds_list)

	#remove high ranked names from list
	TaxID_input_list_str = turn_list_to_DB_list((TaxID_input_list))
	#Set query according to user descendent level:
	#if context.UserFlags_dict['DescendantLevel'] != 'specieLevel':
	query = "select * from ncbi_names_db where TaxID IN (%s) and Type is 'scientific name' and rank IN ('species','subspecies','varietas')" %TaxID_input_list_str
	#else:
	#	query = "select * from ncbi_names_db where TaxID IN (%s) and Type is 'scientific name' and rank IN ('species')" %TaxID_input_list_str
	db_rows = query_ncbi_db(query,db_name)
	TaxId_list,Names_list,TaxId_rank_dict,Name_type_dict,TaxID_Name_dict,TaxId_Parent_dict = extract_names_from_rows(context,db_rows,db_name)

	logger.debug("No High rank in list")
	logger.debug(TaxId_list)
	#sys.exit()

	#get sientific names for TaxIds:
	#keys_list=[]
	#for key in TaxId_rank_dict.keys():
	#	keys_list.append(key)
	#keys_list_for=turn_list_to_DB_list(keys_list)
	tax_scient_str = turn_list_to_DB_list(TaxId_list)
	query = "select * from ncbi_names_db where TaxID IN (%s) and Type is 'scientific name'" %tax_scient_str
	db_rows = query_ncbi_db(query,db_name)
	TaxId_list,Names_list,TaxId_rank_dict,Name_type_dict,TaxID_Name_dict,TaxId_Parent_dict = extract_names_from_rows(context,db_rows,db_name)

	#For debug Only:
	with open(context.working_dir+'/outNR.csv','w') as f_out:
		for key in TaxId_rank_dict.keys():
			#f_out.write('%s,%s,%s\n' %(str(key),str(TaxID_Name_dict[key]).encode("latin-1"),str(TaxId_rank_dict[key])))
			f_out.write('%s,%s,%s\n' %(str(key),str(TaxID_Name_dict[key]),str(TaxId_rank_dict[key])))
			context.species_names_by_tax_id[str(key)] = str(TaxID_Name_dict[key])
			writer.writerow([str(key), str(TaxID_Name_dict[key])])
	f_out.close()

	#Set the lists for OTT pipeline:
	context.species_list_ids = TaxId_list[:]
	context.species_list_names = Names_list[:]
	#context.parent_species_tax_id_by_subsp_tax_id
	#context.parent_species_names_by_subsp_tax_id
	logger.debug("NEW NEW NEW:")
	logger.debug(context.species_list_ids)
	logger.debug(context.species_list_names)
	#Once we have the TaxID list we need to extract all subsp and their parents for merging:
	logger.debug("SYNONYM:")
	logger.debug(context.syn_acc_TaxID_dict)
	logger.debug(context.syn_acc_dict)

	return


#--------------------------------------------------def-------------------------------------------------------
def return_species_list(input_list,context):
# Desc: this function will create the species list and all TaxId data according to user params:
# enables NCBI db alone, spelling correction and Acepted names matching

	# for Accepted names match need to prepare the correct dictionaries:
	# 1. for subsp merge:
	#			parent_tax_id = context.parent_species_tax_id_by_subsp_tax_id[original_tax_id]
	#			parent_name = context.parent_species_names_by_subsp_tax_id[original_tax_id]
	# 2. for Acc/Synonym merge:
	#			context.syn_acc_dict[original_tax_id]
	#			context.accepted_species_names_by_synonym_id[original_tax_id]

	#Get level of retrival:
	depth_level = 'species' #un the future can be either: 'species' (default) or 'subsp'

	names_taxid_f=open(context.species_names_taxid_file,mode="w",encoding='utf8',newline='')
	writer = csv.writer(names_taxid_f,delimiter=',')
	writer.writerow(['TaxID', 'Species_Name'])

	# Decide which DataBase to use for Name matching & Name resolution
	if context.UserFlags_dict['NR_DB_name'] == 'NCBI': #default
		db_name = '/groups/itay_mayrose/share/ploidb/Ploidb_CODE/shared_files/NCBI_Names/NCBI_Names_for_NR_July2017.db'
	elif context.UserFlags_dict['NR_DB_name'] == 'PLT': # the plant list
		db_name = ploidb_config['general']['TPL_ALL_DB']
	else:
		logger.debug("Wrong DataBase for Name resolution")
		sys.exit()
		return


	TaxID_Name_dict=dict()
	TaxID_input_list=[]
	names_for_nr=[]
	Names_input_list=[]
	Names_AfterNR_list=[]
	Request_TaxIds_list=[]

	# Check which inputs are taxIds and which Names
	for item in input_list:
		if item.isdigit():
			TaxID_input_list.append(int(item))
		else:
			Names_input_list.append(item)


	#NAME RESOLUTION section
	#----------------------------------------------------------
	#Handle the names according to user request, only relevant for species names (more than 1 str):
	if context.UserFlags_dict['NameResType'] == 'SpellingCorrection' or context.UserFlags_dict['NameResType'] == 'MatchedNames':
		for name in Names_input_list:
			if ' ' in name:
				names_for_nr.append(name)
		Names_input_list = list(set(Names_input_list) - set(names_for_nr))
		Names_AfterNR_list = db_name_matching(db_name,names_for_nr,context)
		print(Names_AfterNR_list)
	#if context.UserFlags_dict['NameResType'] == 'MatchedNames':
		# Also prepare Syn/Accepted dictionaries:
		#Check_for_Synonyms()
		Names_input_list.extend(Names_AfterNR_list)
	#sys.exit()

	# Convert Names to TaxIDs:
	NamesSTRING_for_query = turn_list_to_DB_list(Names_input_list)
	#query_for_Names = "select * from ncbi_names_db where Name IN (%s) and Type is 'scientific name'" %NamesSTRING_for_query
	query_for_Names = "select * from ncbi_names_db where lower(Name) IN (%s) and Type not in ('common name','authority')" %NamesSTRING_for_query
	db_rows = query_ncbi_db(query_for_Names,db_name)
	TaxId_list,Names_list,TaxId_rank_dict,Name_type_dict,TaxID_Name_dict,TaxId_Parent_dict = extract_names_from_rows(context,db_rows,db_name)
	if TaxID_input_list:
		Request_TaxIds_list = set(TaxID_input_list + TaxId_list)
	else:
		Request_TaxIds_list=list(TaxId_list)
	logger.debug("User requested TaxIds list:")
	logger.debug(Request_TaxIds_list)

	#combine Name HR with TaxId HR
	#Working with Higher ranks
	#by distinguishing between parents with spc kids and parents with higher ranked kids:
	TaxID_input_list=[]
	while Request_TaxIds_list:
		temp_HR_list = Request_TaxIds_list
		Spc_rank_TaxId,High_rank_TaxId = Split_low_high_Ranks(context,temp_HR_list,db_name)
		logger.debug("Spc List - Spc_rank_TaxId:")
		logger.debug(Spc_rank_TaxId)
		logger.debug("Spc List - TaxID_input_list:")
		logger.debug(TaxID_input_list)
		logger.debug("High rank taxa, need to check for kids:")
		logger.debug(High_rank_TaxId)
		#Update theTaxId species list:
		TaxID_input_list.extend(list(set(Spc_rank_TaxId)))
		#Continue with parents taxIds:
		Request_TaxIds_list = list(High_rank_TaxId)

	logger.debug("TaxId list - all")
	logger.debug(Request_TaxIds_list)

	#remove high ranked names from list
	TaxID_input_list_str = turn_list_to_DB_list((TaxID_input_list))
	#Set query according to user descendent level:
	#if context.UserFlags_dict['DescendantLevel'] != 'specieLevel':
	query = "select * from ncbi_names_db where TaxID IN (%s) and Type is 'scientific name' and rank IN ('species','subspecies','varietas')" %TaxID_input_list_str
	#else:
	#	query = "select * from ncbi_names_db where TaxID IN (%s) and Type is 'scientific name' and rank IN ('species')" %TaxID_input_list_str
	db_rows = query_ncbi_db(query,db_name)
	TaxId_list,Names_list,TaxId_rank_dict,Name_type_dict,TaxID_Name_dict,TaxId_Parent_dict = extract_names_from_rows(context,db_rows,db_name)

	logger.debug("No High rank in list")
	logger.debug(TaxId_list)
	#sys.exit()

	#get sientific names for TaxIds:
	#keys_list=[]
	#for key in TaxId_rank_dict.keys():
	#	keys_list.append(key)
	#keys_list_for=turn_list_to_DB_list(keys_list)
	tax_scient_str = turn_list_to_DB_list(TaxId_list)
	query = "select * from ncbi_names_db where TaxID IN (%s) and Type is 'scientific name'" %tax_scient_str
	db_rows = query_ncbi_db(query,db_name)
	TaxId_list,Names_list,TaxId_rank_dict,Name_type_dict,TaxID_Name_dict,TaxId_Parent_dict = extract_names_from_rows(context,db_rows,db_name)

	#For debug Only:
	with open(context.working_dir+'/outNR.csv','w') as f_out:
		for key in TaxId_rank_dict.keys():
			#f_out.write('%s,%s,%s\n' %(str(key),str(TaxID_Name_dict[key]).encode("latin-1"),str(TaxId_rank_dict[key])))
			f_out.write('%s,%s,%s\n' %(str(key),str(TaxID_Name_dict[key]),str(TaxId_rank_dict[key])))
			context.species_names_by_tax_id[str(key)] = str(TaxID_Name_dict[key])
			writer.writerow([str(key), str(TaxID_Name_dict[key])])
	f_out.close()

	#Set the lists for OTT pipeline:
	context.species_list_ids = TaxId_list[:]
	context.species_list_names = Names_list[:]
	#context.parent_species_tax_id_by_subsp_tax_id
	#context.parent_species_names_by_subsp_tax_id
	logger.debug("NEW NEW NEW:")
	logger.debug(context.species_list_ids)
	logger.debug(context.species_list_names)
	#Once we have the TaxID list we need to extract all subsp and their parents for merging:
	logger.debug("SYNONYM:")
	logger.debug(context.syn_acc_TaxID_dict)
	logger.debug(context.syn_acc_dict)

	return