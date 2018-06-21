import operator
from ott_objects_defs.ploidbCommon import *
from ott_objects_defs.ploidbUtils import *
from ott_objects_defs.OutGroupSelectionContext import OutGroupSelectionContext


# Because of long lines in clusters files and accession tables we enlarged the size to enable csv reader to read all lines:
csv.field_size_limit(2000000)

"""
This file contains several methods dealing with finding appropriate outgroup either for a single cluster (locus) or for multiple clusters (loci)
Methods here are called from buildTaxaTree.py
"""

__author__ = 'moshe'

blast_results_header = "outgroup_seq_id,taxonomy_last,organism,gi,gb,sacc,saccver,evalue,bitscore,score,pident,ppos,qlen,slen,length"

#This function will run the check for selecting clusters to be included in the final concatenated alignment:
def select_clusters_for_final_concat(ploidb_context,outgroup_selection):
	# ------------------------------------ Additional Code for Cluster Selection (michal) ---------------------------------------#
	# New Addition of user putgroup -> need to check if included in one of the selected clusters, otherwise notify the user:

	max_concat_seq_length = int(ott_config['common outgroup']['max_concat_seq_length'])
	logger.debug("Maximum concat length is %s" % str(max_concat_seq_length))

	# Initialize list and params for concat clusters:
	Chosen_cluster_ids_list = []
	max_cluster_index_Id_dict = {}
	all_cluster_length = 0
	Chosen_cluster_list = []
	ref_species_list = []
	num_of_clusters = len(ploidb_context.cluster_contexts)
	i = 0
	# This section will determine which clusters will be added to the final concat:
	# First it will add clusters that add the most data (according to wheight as defined below)
	# Then it will continue to add clusters with wheight 0 according to the NumOf taxa in them (max first)
	# it will stop when the sum of median length of clusters is less then 20000
	while all_cluster_length < max_concat_seq_length:  # or i < num_of_clusters:
		added_species_dict = {}
		weight_dict = {}
		length_dict = {}
		for cluster_context in ploidb_context.cluster_contexts:
			max_cluster_index_Id_dict[cluster_context.index] = cluster_context.cluster_id
			if cluster_context.index in Chosen_cluster_list:
				continue
			else:
				cluster_len = get_cluster_ids_total_seq_len(ploidb_context, [cluster_context.cluster_id])
				current_context = ploidb_context.get_cluster_by_id(cluster_context.cluster_id)
				logger.debug(cluster_context.index)
				# NEED to add filter for cluster selection in case of no outgroup:
				# def get_cluster_contribution(outgroup_context, ploidb_context, context_to_check):
				cluster_species_list = current_context.get_list_of_species()
				logger.debug("Checking if to add the following cluster %s, length = %s, species count = %s" % (
				ploidb_context.get_cluster_by_id(cluster_context.cluster_id), str(cluster_len),
				str(len(cluster_species_list))))
				logger.debug("Chosen Cluster List: %s" % Chosen_cluster_list)
				logger.debug("Species List: %s" % cluster_species_list)
				added_species = list(set(cluster_species_list) - set(ref_species_list))
				logger.debug("diff_species List: %s" % added_species)
				weight_dict[cluster_context.index] = len(cluster_species_list) * len(
					added_species)  # Weight is the number of species in cluster TIMES the number os New species it adds to the concat
				logger.debug("Weight %s" % str(weight_dict[cluster_context.index]))
				length_dict[cluster_context.index] = cluster_len
				added_species_dict[cluster_context.index] = added_species
				logger.debug(weight_dict)
				logger.debug(length_dict)
		i += 1
		max_cluster_index = max(weight_dict, key=weight_dict.get)
		all_cluster_length += length_dict[max_cluster_index]
		Chosen_cluster_list.append(max_cluster_index)
		Chosen_cluster_ids_list.append(max_cluster_index_Id_dict[max_cluster_index])
		ref_species_list.extend(added_species_dict[max_cluster_index])
		logger.debug("Updated clusters list: %s" % Chosen_cluster_list)
		logger.debug("Added species list: %s" % added_species_dict[max_cluster_index])
		logger.debug("Adding Cluster %s, Total length: %s" % (str(max_cluster_index), str(all_cluster_length)))
		# Check user outgroup:
		if 'Outgroup_User' in ploidb_context.UserFlags_dict:
			for taxID_str in added_species_dict[max_cluster_index]:
				if taxID_str == ploidb_context.UserOutgroupTaxId:
					logger.debug("Found User Outgroup in added cluster !!!")
					ploidb_context.UserOutgroupInc = 'YES'
				# if str(ploidb_context.UserOutgroupTaxId) in added_species_dict[max_cluster_index]:
				#	logger.debug("Found User Outgroup in added cluster !!!")
				#	ploidb_context.UserOutgroupInc = 'YES'
		if len(Chosen_cluster_list) == num_of_clusters:
			break

	# Update summary file in case User outgroup is not included in the output:
	if 'Outgroup_User' in ploidb_context.UserFlags_dict and ploidb_context.UserOutgroupInc == 'NO':
		with open(ploidb_context.summary_file, 'w') as sum_f:
			sum_f.write("User Outgroup (%s) is Not included in the final Alignment\n" % ploidb_context.UserFlags_dict[
				'Outgroup_User'])
			ploidb_context.UserOutgroupInMSA = False
			#statusFail_LogFile(ploidb_context, 'The User Outgroup you specified was not included in the final MSA and so OneTwoTree did not reconstruct a phylogeny.')
			#raise Exception("Failed - The User Outgroup you specified was not included in the final MSA and so OneTwoTree did not reconstruct a phylogeny.")
	elif 'Outgroup_User' in ploidb_context.UserFlags_dict:
		# Write user outgroup name to file
		with open(ploidb_context.outgroup_file, 'w') as out_f:
			out_f.write(ploidb_context.UserFlags_dict['Outgroup_User'])
		ploidb_context.outgroupSelection = get_selected_outgroup(ploidb_context)
		with open(ploidb_context.summary_file, 'a') as f_sum:
			f_sum.write("User Selected Outgroup: %s\n" % ploidb_context.outgroupSelection)
	elif outgroup_selection == 'None':
		with open(ploidb_context.outgroup_file, 'w') as out_f:
			out_f.write('None')
	# Copy chosen clusters to final list:
	logger.debug("Chosen Cluster index: %s" % Chosen_cluster_list)
	logger.debug("Chosen Cluster Ids: %s " % Chosen_cluster_ids_list)
	ploidb_context.cluster_contexts_for_concat_tree = ploidb_context.get_clusters_by_ids(Chosen_cluster_ids_list)
	logger.debug(ploidb_context.cluster_contexts_for_concat_tree)

	# for cluster_context in ploidb_context.cluster_contexts:
	for cluster_context in ploidb_context.cluster_contexts_for_concat_tree:
		if cluster_context.index in Chosen_cluster_list:
			logger.debug("Adding Cluster %s to concat list" % str(cluster_context.index))
			cluster_context.is_used_for_concat_tree = True

	return# Chosen_cluster_list

def init_cluster_db(cluster_context, context,outgroup_selection):
	logger.debug("length of species_list - " + str(len(context.species_list_names)))
	logger.debug("species list: " + ', '.join(context.species_list_names))
	if not cluster_context.is_db_init:
		cache_flag = ott_config.getboolean('general', 'use_cache')
		# raise Exception("killing for outgroup in init_cluster_db")
		logger.info(
			"Writing outgroup for the cluster in  %s" % cluster_context.all_seqs_fasta_filename_no_multiple_accessions)
		#Skip this blast in case of No Outgroup - not needed, saving running time (michal):
		if context.UserFlags_dict['Outgroup_Flag'] == 'Single':
			logger.debug("User selected a Single Outgroup, perform blastn per cluster for outgroup retrieval ")
			blast_sequences_into_db_local(cluster_context, context.db_file,
								cluster_context.all_seqs_fasta_filename_no_multiple_accessions,
								cluster_context.cluter_work_dir, cache_flag, context,outgroup_selection)
		else:
			#Generate the seq-to-blast (representative seq) to be used with Large taxa addition:
			clusterRepresentativeSequence = cluster_context.get_cluster_seq_to_blast()
			seqToBlastFilename = cluster_context.cluter_work_dir + "/seq-to-blast"
			SeqIO.write(clusterRepresentativeSequence, seqToBlastFilename, "fasta")
			logger.debug("User selected none or user Outgroup, skipping blastn per cluster") #  skipping - blast_sequences_into_db

		cluster_context.locus_tree_outgroup_selection.blast_results_outside_cluster = number_of_blast_results_outside_cluster(
			context.db_file, cluster_context.cluster_id)
		min_blast_results_outside_genus = int(ott_config['common outgroup']['min_blast_results_outside_genus'])
		if cluster_context.locus_tree_outgroup_selection.blast_results_outside_cluster < min_blast_results_outside_genus:
			logger.info(
				"BLAST returned no results that can be used as outgroup for cluster %s" % cluster_context.cluster_id)
		else:
			logger.info(
				"BLAST returned %d results that are not within the genus" % cluster_context.locus_tree_outgroup_selection.blast_results_outside_cluster)

			try:
				out_seqs = write_single_locus_outgroup_to_file(context, cluster_context, context.id)
				if out_seqs is None or len(out_seqs) == 0:
					logger.warning(
						"Failed to find single locus outgroup for cluster %s. Cluster will be skipped" % cluster_context.cluster_id)
				else:
					cluster_context.locus_tree_outgroup_selection.outgroup_found = True
			except Exception as e:
				logger.exception(
					"Exception caught while trying to set outgroup to cluster %s" % cluster_context.cluster_id)

		cluster_context.is_db_init = True
	else:
		logger.info("cluster context %s blast db already initialized. Skipping" % cluster_context.cluster_id)



#This function is used in the new method of using blast on a local database:
# The database is created from the NCBI genbank db file
def prepare_blast_results_for_insert_db(blastResultFilename, blast_results_for_insert_filename, outgroup_seq_id, context):
	logger.debug("Preparing blast results for insert")
	logger.debug("##################################")
	#gbDict = context.get_gb_seqs_dict()

	with open(blastResultFilename, 'r') as in_file, open(blast_results_for_insert_filename, 'w') as out_file:
		writer = csv.writer(out_file, delimiter=',', dialect=csv.excel)
		writer.writerow(blast_results_header.split(","))
		# out_file.write(blast_results_header + "\n")
		missing_saccver = 0
		for line in in_file:
			line_words = line.rstrip().split(',')  #rstrip() removes "\n" so last column won't be considered a string
			saccver = line_words[1].split('.')[0]
			sseqid = line_words[0]
			sseqid_words = sseqid.split("|")
#			if saccver not in gbDict:
#				missing_saccver += 1
#				logger.debug("Couldn't find key %s in seq dictionary. It will be ignored in the BLAST results", saccver)
#			else:
			gi = sseqid_words[0]  # Taxaon ID
			gb = sseqid_words[1]  # Accession number
			context.outgroup_AccessionkeySeq_dict[gb] = sseqid
			sacc = gb.split('.')[0]
			saccver=gb
			organism = sseqid_words[2].replace('_',' ')
			genus = organism.split("_")[0]
			if not " " in organism:
				logger.warning(
					"Organism %s (gi=%s) doesn't contain at least 2 words and will be skipped from blast results" % (
						organism, gi))
			else:
				writer.writerow([outgroup_seq_id, genus, organism, gi, gb, sacc, saccver, line_words[3],
								 line_words[4], line_words[5], line_words[6], line_words[7], line_words[8],
								 line_words[9], line_words[10]])

		if missing_saccver > 0:
			logger.warn("%d keys were not found in the dictionary" % missing_saccver)

	logger.debug("DONE preparing blast results for insert")


def insert_blast_results_to_db(db_file, outgroup_seq_id, csv_file):
	logger.debug("Inserting BLAST results in the DB")
	conn = get_db_conn(db_file)
	curs = conn.cursor()

	with open(csv_file, 'r') as infile:
		# csv.DictReader uses first line in file for column headings by default
		dr = csv.DictReader(infile, delimiter=',', dialect=csv.excel)
		to_db = [(i['outgroup_seq_id'], i['taxonomy_last'], i['organism'], i['gi'], i['gb'], i['sacc'], i['saccver'],
				  i['evalue'], i['bitscore'], i['score'], i['pident'], i['ppos'], i['qlen'], i['slen'], i['length']) for
				 i in dr]
	curs.executemany("INSERT INTO blast_results (" + blast_results_header + ") VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);",
					 to_db)
	# Save (commit) the changes
	conn.commit()
	# We can also close the connection if we are done with it.
	# Just be sure any changes have been committed or they will be lost.
	conn.close()

def write_filtered_results_to_csv(db_file, csv_filename):
	with get_db_conn(db_file) as conn:
		cursor = conn.cursor()
		cursor.execute("select * from blast_results_filtered;")

		csv_writer = csv.writer(open(csv_filename, "wt"))
		csv_writer.writerow([i[0] for i in cursor.description])  # write headers
		csv_writer.writerows(cursor)
		del csv_writer  # this will close the CSV file


# Yoni - changing the outgroup selection:
#        instead of setting is_seq_in_species_list = 1 when organism contains genus name
#        checking is organism is in context.species_list
#        is_seq_in_species_list column used to be is_seq_in_genus
def create_view_for_cluster(db_file, cluster_context, outgroup_seq_id, fasta_filename, species_list):

	logger.debug("Creating view for cluster %s" % cluster_context.cluster_id)
	with get_db_conn(db_file) as conn:
		# Getting the max evalue of all of the seqs
		curs = conn.cursor()

		out_query = "SELECT max(evalue) as eval FROM blast_results_filtered where outgroup_seq_id = '%s'" % outgroup_seq_id
		logger.debug("Getting max evalue without restricting GI using the following query %s ", out_query)
		logger.debug("Using Query: SELECT max(evalue) as eval FROM blast_results_filtered where outgroup_seq_id = '%s'" % outgroup_seq_id)
		curs.execute(out_query)
		rows = curs.fetchall()

		for row in rows:
			eval = row[0]
			logger.debug("max evalue for all is %s" % str(eval))

		# Getting the max evalue without the genus GIs
		curs = conn.cursor()

		gi_list_str = gi_list_str_from_fasta_file(fasta_filename)
		logger.debug("Print list of accession numbers (as gi's) gi_list_str")
		logger.debug(gi_list_str)
		view_name = cluster_context.get_view_name()
		try:
			curs.execute("DROP VIEW %s" % view_name)
		except sqlite3.OperationalError:
			# View doesn't exists
			pass
		update_cluster_gi_sql = "UPDATE blast_results SET is_cluster_gi = 1 where outgroup_seq_id = '%s' AND gi in (%s)" % (
			outgroup_seq_id, gi_list_str)
		curs.execute(update_cluster_gi_sql)
		logger.debug("Update Database: UPDATE blast_results SET is_cluster_gi = 1 where outgroup_seq_id = %s AND gi in (%s)" % (outgroup_seq_id, gi_list_str))

		# yoni - creating new update sql command
		logger.debug("length of species_list - " + str(len(species_list)))
		logger.debug("species list: " + ', '.join(species_list))

		# If species list is too long (> 1000) than needto split the sql update:
		if len(species_list) > 800:
			logger.debug("Species List too long for sql update, perform update by chanks:")
			chank_species_lists=[species_list[x:x+800] for x in range(0, len(species_list), 800)]
			for chank_species_list in chank_species_lists:
				logger.debug("Chank species list: " + ', '.join(chank_species_list))
				placeholder = '?'
				placeholders = ', '.join(placeholder * len(chank_species_list))
				update_is_seq_in_species_list_sql = "UPDATE blast_results SET is_seq_in_species_list = 1 where outgroup_seq_id = '%s' and organism IN (%s)" % (outgroup_seq_id, placeholders)
				logger.debug("Update sql command is- %s" % update_is_seq_in_species_list_sql)
				curs.execute(update_is_seq_in_species_list_sql, chank_species_list)
		else:
			placeholder = '?'
			placeholders = ', '.join(placeholder * len(species_list))
			update_is_seq_in_species_list_sql = "UPDATE blast_results SET is_seq_in_species_list = 1 where outgroup_seq_id = '%s' and organism IN (%s)" % (outgroup_seq_id, placeholders)
			logger.debug("Update sql command is- %s" % update_is_seq_in_species_list_sql)
			curs.execute(update_is_seq_in_species_list_sql, species_list)

		update_is_min_score_for_acc_sql = "update blast_results " \
										  "set is_min_score_for_acc = '1' " \
										  " where  (gi||'#' ||cast(bitscore as string)) in (select gi||'#' ||cast(bitscore as string) from blast_results_max);"
		logger.debug(update_is_min_score_for_acc_sql)
		logger.debug("Print for Debug: update_is_min_score_for_acc_sql")
		curs.execute(update_is_min_score_for_acc_sql)
		evalue_threshold = float(ott_config['common outgroup']['evalue_threshold'])
		create_view_sql = "CREATE VIEW %s AS select * from blast_results_filtered where outgroup_seq_id = '%s' AND is_cluster_gi = 1 AND evalue > %f" % (
			view_name, outgroup_seq_id, evalue_threshold)
		logger.debug("Print for Debug: create_view_sql")
		logger.debug(create_view_sql)
		curs.execute(create_view_sql)
		logger.debug("Created the following view %s " % create_view_sql)

		verify_outgroup_exists = " select count(*) from blast_results_filtered where outgroup_seq_id = '%s' AND is_cluster_gi = 0 AND is_seq_in_species_list = 0" % outgroup_seq_id
		logger.debug("Verifying BLAST results contains outgroup using %s " % verify_outgroup_exists)
		curs.execute(verify_outgroup_exists)
		result = curs.fetchone()
		number_of_rows = result[0]
		logger.debug("BLAST returned %d results that are not within the genus" % number_of_rows)

"""
Checking if BLAST returned at least one sequence that is outside the cluster, and not part of this genus
This is just a naiive check since when selecting outgroup, there are more restrictions
"""

def number_of_blast_results_outside_cluster(db_file, outgroup_seq_id):
	with get_db_conn(db_file) as conn:
		# Getting the max evalue of all of the seqs
		curs = conn.cursor()

		verify_outgroup_exists = " select count(*) from blast_results_filtered where outgroup_seq_id = '%s' AND is_cluster_gi = 0 AND is_seq_in_species_list = 0" % outgroup_seq_id
		logger.debug("Verifying BLAST results contains outgroup using %s " % verify_outgroup_exists)
		curs.execute(verify_outgroup_exists)
		result = curs.fetchone()
		number_of_rows = result[0]
		return number_of_rows


"""
Calculate mean and standard deviation of data x[]:
	mean = {\sum_i x_i \over n}
	std = sqrt(\sum_i (x_i - mean)^2 \over n-1)
"""


def meanstdv(x):
	from math import sqrt

	n, mean, std = len(x), 0, 0
	for a in x:
		mean = mean + a
	mean /= float(n)
	for a in x:
		std += (a - mean) ** 2
	std = sqrt(std / float(n - 1))
	return std


def init_cluster_blast_stats(db_file, cluster_context, genus_name):
	"""
	@type cluster_context: ClusterContext
	"""
	entire_cluster_where_clause = "outgroup_seq_id = '%s'" % cluster_context.cluster_id
	init_cluster_blast_stats_per_attr(db_file=db_file,
									  stat_summary_data=cluster_context.locus_tree_outgroup_selection.cluster_stats,
									  results_sql_clause=entire_cluster_where_clause, colname="evalue")
	init_cluster_blast_stats_per_attr(db_file=db_file,
									  stat_summary_data=cluster_context.locus_tree_outgroup_selection.cluster_stats,
									  results_sql_clause=entire_cluster_where_clause, colname="bitscore")

	genus_cluster_where_clause = "outgroup_seq_id = '%s' and is_seq_in_species_list = 1" % cluster_context.cluster_id
	init_cluster_blast_stats_per_attr(db_file=db_file,
									  stat_summary_data=cluster_context.locus_tree_outgroup_selection.genus_stats,
									  results_sql_clause=genus_cluster_where_clause, colname="evalue")
	init_cluster_blast_stats_per_attr(db_file=db_file,
									  stat_summary_data=cluster_context.locus_tree_outgroup_selection.genus_stats,
									  results_sql_clause=genus_cluster_where_clause, colname="bitscore")

	non_genus_cluster_where_clause = "outgroup_seq_id = '%s' and is_seq_in_species_list = 0" % cluster_context.cluster_id
	init_cluster_blast_stats_per_attr(db_file=db_file,
									  stat_summary_data=cluster_context.locus_tree_outgroup_selection.non_genus_stats,
									  results_sql_clause=non_genus_cluster_where_clause, colname="evalue")
	init_cluster_blast_stats_per_attr(db_file=db_file,
									  stat_summary_data=cluster_context.locus_tree_outgroup_selection.non_genus_stats,
									  results_sql_clause=non_genus_cluster_where_clause, colname="bitscore")


def init_cluster_blast_stats_per_attr(db_file, stat_summary_data, results_sql_clause, colname):
	"""
	@type cluster_context: ClusterContext
	@type stat_summary_data: BlastStatsSummary
	"""
	with get_db_conn(db_file) as conn:
		curs = conn.cursor()

		out_query = "SELECT max(%s) as maxval,min(%s) as minval,avg(%s) as avgval " \
					"FROM blast_results_filtered " \
					"WHERE %s " % (colname, colname, colname, results_sql_clause)
		logger.debug("Getting value stats from %s for column %s using the following query %s " % (db_file, colname, out_query))
		curs.execute(out_query)
		rows = curs.fetchall()

		if len(rows) != 1:
			logger.critical("Expected single row and got %d rows" % len(rows))
		else:
			row = rows[0]
			try:
				if colname == "evalue":
					stat_summary_data.max_eval = float(row[0])
					stat_summary_data.min_eval = float(row[1])
					stat_summary_data.avg_eval = float(row[2])
				elif colname == "bitscore":
					stat_summary_data.max_score = float(row[0])
					stat_summary_data.min_score = float(row[1])
					stat_summary_data.avg_score = float(row[2])
				else:
					raise Exception("Unknown column name %s" % colname)
			except(TypeError, ValueError):
				logger.error("Failed to find outgroup using %s. min/max didn't return float value " % colname)

		out_query = "SELECT %s" \
					" FROM blast_results_filtered " \
					" WHERE %s " \
					" ORDER BY %s ASC " % (colname, results_sql_clause, colname)
		logger.debug("Getting %s value stddev using the following query %s " % (colname, out_query))
		curs.execute(out_query)
		rows = curs.fetchall()
		stat_summary_data.number_of_results = len(rows)
		if stat_summary_data.number_of_results == 0:
			logger.debug("Zero results for blast stats of %s. Returning..." % (stat_summary_data.desc))

		logger.debug("There are %d results for blast stats of %s " % (
			stat_summary_data.number_of_results, stat_summary_data.desc))
		all_vals = list()
		for row in rows:
			all_vals.append(float(row[0]))
		percentiles = list()
		for percent in range(0, 100, 1):
			index_of_percentile = int(len(all_vals) * percent / 100)
			# logger.debug("Adding to index %d the value %d" % (index_of_percentile,all_vals[index_of_percentile]))
			percentiles.append(all_vals[index_of_percentile])

		low_precentile = float(ott_config['common outgroup']['low_precentile'])
		index_of_low_percentile_val = int(len(all_vals) * low_precentile)
		index_of_high_percentile_val = int(len(all_vals) * (1 - low_precentile))
		if len(all_vals) <= 1:
			stddev = 0
		else:
			stddev = meanstdv(all_vals)

		logger.debug("stdev is %f. low and high indexes are [%d,%d] for blast stats of %s " % (stddev,
																							   index_of_low_percentile_val,
																							   index_of_high_percentile_val,
																							   stat_summary_data.desc))

		if colname == "evalue":
			stat_summary_data.stddev_eval = stddev
			stat_summary_data.low_percentile_eval = all_vals[index_of_low_percentile_val]
			stat_summary_data.high_percentile_eval = all_vals[index_of_high_percentile_val]
			stat_summary_data.percentile_eval = percentiles
			stat_summary_data.percentile_eval_str = ",".join(str(x) for x in percentiles)
		elif colname == "bitscore":
			stat_summary_data.stddev_score = stddev
			stat_summary_data.low_percentile_score = all_vals[index_of_low_percentile_val]
			stat_summary_data.high_percentile_score = all_vals[index_of_high_percentile_val]
			stat_summary_data.percentile_score = percentiles
			stat_summary_data.percentile_score_str = ",".join(str(x) for x in percentiles)
		else:
			raise Exception("Unknown column name %s" % colname)


def init_cluster_locus_outgroup_params(cluster_context):
	"""
	@type cluster_context: ClusterContext
	"""
	genus_min_score = cluster_context.locus_tree_outgroup_selection.genus_stats.min_score
	genus_low_percentile_score = cluster_context.locus_tree_outgroup_selection.genus_stats.percentile_score[5]

	non_genus_avg_score = cluster_context.locus_tree_outgroup_selection.non_genus_stats.avg_score
	non_genus_high_percentile_score = cluster_context.locus_tree_outgroup_selection.non_genus_stats.percentile_score[90]
	non_genus_max_score = cluster_context.locus_tree_outgroup_selection.non_genus_stats.max_score
	logger.debug("Init cluster locus outgroup parameters. [genus: min_score = %f, low_percentile=%f] "
				 "[non_genus: avg_score = %f max_score = %f high_percentile = %f ]" %
				 (genus_min_score,
				  genus_low_percentile_score,
				  non_genus_avg_score,
				  non_genus_max_score,
				  non_genus_high_percentile_score))

	genus_max_eval = cluster_context.locus_tree_outgroup_selection.genus_stats.max_eval

	cluster_context.locus_tree_outgroup_selection.minscore = genus_min_score
	# cluster_context.locus_tree_outgroup_selection.eval = genus_max_eval
	if non_genus_max_score <= genus_low_percentile_score:
		logger.info("SUCCESS FINDING OUTGROUP :) non_genus_max_score <= genus_min_score")
	elif non_genus_high_percentile_score < genus_low_percentile_score:
		cluster_context.locus_tree_outgroup_selection.minscore = genus_low_percentile_score
		logger.info(
			"partial SUCCESS FINDING OUTGROUP non_genus_high_percentile_score < genus_low_percentile_score [%f, %f] Using genus_low_percentile_score" % (
				non_genus_high_percentile_score, genus_low_percentile_score))
	else:
		if genus_min_score < non_genus_avg_score:
			logger.warning("Genus min score is lower than average score of non genus scores")
			# In this case this locus should not be used
			return False
		else:  # genus_min_score > non_genus_avg_score
			cluster_context.locus_tree_outgroup_selection.minscore = genus_low_percentile_score
			logger.warning("Genus score is mixed up with non genus scores.Using genus_low_percentile_score")

	return True

# if genus_max_eval > EVALUE_THRESHOLD:
#		logger.error("Evalue that will be used for outgroup is %f which is higher than 10e-50. Outgroup is not optimal " % genus_max_eval)


#	logger.debug("Init eval params. genus_max_eval = %f and non_genus_min_eval = %f" % (genus_max_eval,non_genus_min_eval))
#	if genus_max_eval <= non_genus_min_eval :
#		logger.debug("genus_max_eval <= non_genus_min_eval. Setting eval to non_genus_min_eval ")
#		cluster_context.locus_tree_outgroup_selection.eval= non_genus_min_eval
#	else:
#		logger.warning("max eval of genus (%f) intersects with non genus eval(%f). Using low percentile eval" % (genus_max_eval, non_genus_min_eval))
#		cluster_context.locus_tree_outgroup_selection.eval= cluster_context.locus_tree_outgroup_selection.genus_stats.high_percentile_eval


"""
The purpose
"""

def get_seqs_by_key(context,keys, new_organism=None):
	#if not ploidb_context.gbDict:   # made this change to save running time (to create it just once)
	#gbDict = get_indexed_seq_hash(ott_config['general']['genbankDir'])
	if context.gb_seqs_dict is not None:
		gbDict = context.get_gb_seqs_dict()
	else:
		gbDict = context.gb_seqs_dict()
	out_seqs = list()
	for key in keys:
		out_seq_saccver = str(key)
		logger.debug("saccver key is %s" % out_seq_saccver)
		if out_seq_saccver in gbDict:
			seq = gbDict[out_seq_saccver]

			standarizeSequence(seq=seq, override_organism=new_organism)
			# TODO: use more than a single outgroup
			if len(out_seqs) < 1:
				logger.debug("adding the following seq as out group " + out_seq_saccver)
				out_seqs.append(seq)
			else:
				return out_seqs
	return out_seqs

#This functin will return the sequnce fom the blast database according to its key:
def get_seqs_by_key_db(context,keys, new_organism=None):

	blastDbPath = ott_config['general']['DB_DIR'] + ott_config['concat_headir']['blastDbPath']
	out_seqs = list()
	for key in keys:
		seq_identifier = context.outgroup_AccessionkeySeq_dict[key]
		logger.debug("blastdbcmd -> -db %s -entry %s" %(blastDbPath,seq_identifier))
		out_data = context.working_dir+'/out_temp.fasta'
		os.system('blastdbcmd -db %s -entry "%s" > %s' %(blastDbPath,seq_identifier,out_data))
		for seq in SeqIO.parse(out_data, "fasta"):
			out_seq_saccver = str(key)
			logger.debug("saccver key is %s" % out_seq_saccver)
			desc= seq.description
			gi = desc.split('|')[0]
			taxon_id = desc.split('|')[1]
			organism = desc.split('|')[2]
			seq.description = "gi|%s|taxonid|%s|organism|%s" % (gi,taxon_id,organism)
			logger.debug("Outgroup description line: %s" %seq.description)
			#standarizeSequence(seq=seq, override_organism=new_organism)
			# TODO: use more than a single outgroup
			if len(out_seqs) < 1:
				logger.debug("adding the following seq as out group " + out_seq_saccver)
				out_seqs.append(seq)
			else:
				return out_seqs
	return out_seqs

def get_common_outgroup_from_db_results(ploidb_context,db_file, outgroup_seq_id, eval, minscore, outgroup_name, outgroup_type):
	with get_db_conn(db_file) as conn:
		curs = conn.cursor()

		extra_conditions = ""
		if eval is not None and minscore is not None:
			extra_conditions += " AND evalue >= %f and bitscore * 1.0 <= %f " % (-1, minscore)
		if outgroup_type == "organism":
			extra_conditions += " AND organism = '%s'" % outgroup_name
		elif outgroup_type == "genus":
			extra_conditions += " AND taxonomy_last = '%s'" % outgroup_name
		else:
			logger.error("Unknown outgroup type %s" % outgroup_type)
		logger.debug("adding the following to outgroup query %s" % extra_conditions)

		out_query = " SELECT saccver " \
					" FROM blast_results_filtered " \
					" WHERE outgroup_seq_id = '%s' %s " \
					" order by evalue asc, bitscore desc, gi" % ( outgroup_seq_id, extra_conditions)
		logger.debug("Getting out group using the following query %s " % out_query)
		curs.execute(out_query)
		rows = curs.fetchall()

		keys = list()
		for row in rows:
			keys.append(row[0])
		new_organism = None
		if outgroup_type == "genus":
			new_organism = outgroup_name
		#return get_seqs_by_key(ploidb_context,keys, new_organism)
		return get_seqs_by_key_db(ploidb_context,keys, new_organism)


def get_locus_outgroup_from_db_results(context,db_file, outgroup_seq_id, eval, minscore):
	with get_db_conn(db_file) as conn:
		curs = conn.cursor()

		extra_conditions = ""
		extra_conditions += " AND evalue * 1.0 >= %f and bitscore * 1.0 <= %f " % (eval, minscore)
		logger.debug("adding the following to outgroup query %s" % extra_conditions)

		out_query = " SELECT saccver " \
					" FROM blast_results_filtered " \
					" WHERE is_cluster_gi = 0 and is_seq_in_species_list = 0 AND outgroup_seq_id = '%s' %s " \
					" order by evalue asc, bitscore desc,gi" % (outgroup_seq_id, extra_conditions)
		logger.debug("Getting out group using the following query %s " % out_query)
		curs.execute(out_query)
		rows = curs.fetchall()

		keys = list()
		for row in rows:
			keys.append(row[0])

		#return get_seqs_by_key(context,keys, None)
		return get_seqs_by_key_db(context,keys, None)


def is_blast_result_already_in_db(db_file, outgroup_seq_id):
	with get_db_conn(db_file) as conn:
		curs = conn.cursor()
		out_query = "SELECT * FROM blast_results_filtered where outgroup_seq_id = '%s'" % outgroup_seq_id
		logger.debug("Checking if there are already blast results in db using the following query %s " % out_query)
		curs.execute(out_query)
		rows = curs.fetchall()
		results_exist = len(rows) > 0
		logger.debug("%i results were found for %s returning %s" % (len(rows), outgroup_seq_id, results_exist))
	return results_exist


def get_filename_for_blast_cache(seq, genus_name):
	gi = getPropertyFromFastaSeqHeader(seq.description, "gi")
	path_for_cache = os.path.join(ott_config['general']['cache_path'], genus_name)
	filename_for_cache = os.path.join(path_for_cache, "bl_gi%s" % gi)
	return filename_for_cache

#This should replace the previous function that performed blast vs the nt database.
# now th database is performed vs a local blast database created from all seqs in the NCBI Genbank db
def blast_sequences_into_db_local(cluster_context, db_file, fasta_filename, temp_results_dir, use_cached_results, context,outgroup_selection):
	blastBinPath = ott_config['general']['blastBinPath']
	blastDbPath = ott_config['general']['DB_DIR'] + ott_config['concat_headir']['blastDbPath']
	logger.debug("Blasting seqs into db for %s" % cluster_context)
	genus_name = context.id
	giListOutFilename = temp_results_dir + "/negative-gis-for-outgroup"
	giListFromFastaFile(fasta_filename, giListOutFilename)
	clusterRepresentativeSequence = cluster_context.get_cluster_seq_to_blast()
	filename_for_blast_cache = get_filename_for_blast_cache(clusterRepresentativeSequence, genus_name)
	seqToBlastFilename = temp_results_dir + "/seq-to-blast"
	blastResultFilename = temp_results_dir + "/blast-results"
	SeqIO.write(clusterRepresentativeSequence, seqToBlastFilename, "fasta")
	# TODO: check if this is what should be used as identifier
	outgroup_seq_id = cluster_context.cluster_id
	logger.debug("Using %s as identifier for this cluter", outgroup_seq_id)
	blast_result_already_in_db = is_blast_result_already_in_db(db_file, outgroup_seq_id)
	# This is the quickest cache - BLAST results already in DB
	#if (not use_cached_results or not blast_result_already_in_db):
	if (not blast_result_already_in_db):
		# at least 500 results
		if outgroup_selection == 'None':
			max_results = 10
		else:
			max_results = max(500, cluster_context.num_of_cluster_seqs * 50)
			max_results = min(15000, max_results)
		blastCommand = blastBinPath + "blastn" + " -query " + seqToBlastFilename + " -db " + blastDbPath + \
					   '  -num_threads 8 -max_target_seqs ' + str(
			max_results) + ' -outfmt "10 sseqid sacc saccver evalue bitscore score pident ppos qlen slen length"'
		# This is the second cache - BLAST results are on a file but still need to be prepared and

		if ((not use_cached_results) or (not os.path.exists(filename_for_blast_cache)) or
		#if ((not os.path.exists(filename_for_blast_cache)) or
		 (not is_non_zero_file(filename_for_blast_cache))):
			logger.debug(
				"Running BLAST use_cached_results=%s blast_result_already_in_db=%s filename_for_blast_cache=%s",
				use_cached_results, blast_result_already_in_db, filename_for_blast_cache)
			# TODO: is this a good estimation?
			#logger.debug("the blast files from cache are of size=%s",os.path.getsize(filename_for_blast_cache))

			exec_external_command_redirect_output(blastCommand, blastResultFilename,
												  temp_results_dir + "/blast-err.log")
			# copying results to cache only if using cache
			if use_cached_results:
				shutil.copy(blastResultFilename, filename_for_blast_cache)
		else:
			logger.debug(
				"Skipped BLAST. will NOT execute the following command %s. use_cached_results=%s blast_result_already_in_db=%s filename_for_blast_cache=%s",
				blastCommand, use_cached_results, blast_result_already_in_db, filename_for_blast_cache)

		blast_results_for_insert_filename = temp_results_dir + "/blast-results-for-insert.csv"

		logger.debug("Executing prepare_blast_results_for_insert_db(%s, %s, %s, %s)" %
					 (blastResultFilename, blast_results_for_insert_filename, outgroup_seq_id, context))

		prepare_blast_results_for_insert_db(blastResultFilename, blast_results_for_insert_filename, outgroup_seq_id,
										 context)

		insert_blast_results_to_db(db_file, outgroup_seq_id, blast_results_for_insert_filename)
		create_view_for_cluster(db_file, cluster_context, cluster_context.cluster_id,
								cluster_context.all_seqs_fasta_filename, context.species_list_names)
		write_filtered_results_to_csv(db_file, temp_results_dir + "/blast-results-filtered.csv")
	else:
		logger.debug(
			"Skipped BLAST and processing of the BLAST data. use_cached_results=%s blast_result_already_in_db=%s",
			use_cached_results, blast_result_already_in_db)


def write_single_locus_outgroup_to_file(context, cluster_context, genus_name):
	output_filename = cluster_context.outgroup_seq_full_filename
	#	fasta_filename = cluster_context.all_seqs_fasta_filename
	#	create_view_for_cluster(context.db_file,cluster_context,cluster_context.cluster_id,fasta_filename,genus_name)

	init_cluster_blast_stats(context.db_file, cluster_context, genus_name)
	are_params_ok = init_cluster_locus_outgroup_params(cluster_context)

	if not are_params_ok:
		logger.debug("At: write_single_locus_outgroup_to_file, LINE -> if not are_params_ok:")
		logger.debug(are_params_ok)
		return None
	else:
		eval = cluster_context.locus_tree_outgroup_selection.eval
		minscore = cluster_context.locus_tree_outgroup_selection.minscore

		out_seqs = get_locus_outgroup_from_db_results(context,context.db_file, cluster_context.cluster_id, eval, minscore)
		if output_filename is not None:
			if len(out_seqs) > 0:
				#out_handle = open(output_filename, 'w')
				logger.debug(
					"Writing %i out group seqs the following to output file %s" % (len(out_seqs), output_filename))
				write_standarize_fasta(output_filename, out_seqs)
			#SeqIO.write(out_seqs, out_handle, "fasta")
			#out_handle.close()
		return out_seqs


def get_clusters_lists_per_outgroup(conn, ploidb_context, outgroup_name):
	curs = conn.cursor()

	# Iterating all possible outgroups
	query = "select cluster_id " \
			"from outgroup_clusters " \
			"where outgroup_name='%s'" % escape_str_for_sql_query(outgroup_name)

	curs.execute(query)
	cluster_ids_rows = curs.fetchall()

	clusters_ids_found = list()
	for cluster_ids_row in cluster_ids_rows:
		clusters_ids_found.append(str(cluster_ids_row[0]))

	logger.debug(
		"Executing query %s returned %d results %s" % (query, len(clusters_ids_found), ",".join(clusters_ids_found)))

	clusters_ids_without_outgroup_seq = list()
	clusters_ids_with_outgroup_seq = list()
	for cluster_id in ploidb_context.get_all_clusters_ids():
		if cluster_id in clusters_ids_found:
			clusters_ids_with_outgroup_seq.append(cluster_id)
		else:
			clusters_ids_without_outgroup_seq.append(cluster_id)

	logger.debug("Outgroup %s: clusters with outgroup %s without %s" % (
		outgroup_name, ",".join(clusters_ids_with_outgroup_seq), ",".join(clusters_ids_without_outgroup_seq)))

	return clusters_ids_with_outgroup_seq, clusters_ids_without_outgroup_seq


def get_cluster_contribution(outgroup_context, ploidb_context, context_to_check):
	species_in_new_cluster = context_to_check.get_list_of_species()
	species_in_all_selected_clusters = outgroup_context.get_species_in_all_selected_clusters()
	#spc that are in the prev selected clusters but are not represented in the candidate cluster
	species_not_in_new_cluster = species_in_all_selected_clusters - species_in_new_cluster
	species_only_in_new_cluster = species_in_new_cluster - species_in_all_selected_clusters


	if len(species_not_in_new_cluster) > 0:
		ratio_of_new_to_existing = len(species_only_in_new_cluster) / len(species_not_in_new_cluster)
	else:
		ratio_of_new_to_existing = 1

	species_in_all_selected_clusters.update(species_in_new_cluster)
	if len(species_in_all_selected_clusters) > 0:
		ratio_of_new_out_of_all = len(species_in_new_cluster) / len(species_in_all_selected_clusters)
	else:
		ratio_of_new_out_of_all = 1

	return ratio_of_new_to_existing, ratio_of_new_out_of_all


def check_and_add_cluster_to_concat(cid_to_check, ploidb_context, outgroup_context):
	logger.debug("Checking the following cluster %s" % ploidb_context.get_cluster_by_id(cid_to_check))
	current_outgroup_len = get_cluster_ids_total_seq_len(ploidb_context,
														 outgroup_context.selected_cids_with_outgroup_seq)
	current_total_concat_len = get_cluster_ids_total_seq_len(ploidb_context, outgroup_context.get_selected_cids())

	cluster_to_check_len = get_cluster_ids_total_seq_len(ploidb_context, [cid_to_check])

	new_outgroup_len = current_outgroup_len
	new_total_concat_len = current_total_concat_len
	if cid_to_check in outgroup_context.all_cids_with_outgroup_seq:
		new_outgroup_len += cluster_to_check_len
	new_total_concat_len += cluster_to_check_len

	outgroup_to_concat_seq_ratio = new_outgroup_len / new_total_concat_len

	# New gaps checking
	context_to_check = ploidb_context.get_cluster_by_id(cid_to_check)

	ratio_of_new_to_existing, ratio_of_new_out_of_all = get_cluster_contribution(outgroup_context, ploidb_context,
																				 context_to_check)

	logger.debug(
		"Got the following: cluster_to_check_len=%d current_outgroup_len=%d current_total_concat_len=%d new_outgroup_len=%d new_total_concat_len=%d ratio_of_new_to_existing=%f ratio_of_new_out_of_all=%f" %
		(cluster_to_check_len, current_outgroup_len, current_total_concat_len, new_outgroup_len, new_total_concat_len,
		 ratio_of_new_to_existing, ratio_of_new_out_of_all))

	max_concat_seq_length = int(ott_config['common outgroup']['max_concat_seq_length'])
	min_species_ratio_for_cluster = float(ott_config['common outgroup']['min_species_ratio_for_cluster'])
	new_to_existing_species_threshold_ratio = float(
		ott_config['common outgroup']['new_to_existing_species_threshold_ratio'])
	min_coverage = float(ott_config['common outgroup']['min_coverage'])
	if max_concat_seq_length < new_total_concat_len:
		logger.debug("Checking outgroup %s cluster_id %s - new_total_concat_len is too high %d returning false" % (
			outgroup_context.outgroup_name, cid_to_check, new_total_concat_len))
		# Reached max length - no way this cluster will be choosen
		outgroup_context.add_cluster_to_not_selected(cid_to_check)
		return False
	if outgroup_to_concat_seq_ratio < min_coverage:
		logger.info(
			"Checking outgroup %s cluster_id %s - outgroup_to_concat_seq_ratio is too low %f returning false" % (
				outgroup_context.outgroup_name, cid_to_check, outgroup_to_concat_seq_ratio))
		if len(outgroup_context.candidate_cids_with_outgroup_seq) > 0:
			logger.debug(
				"Since there are still candidate clusters with outgroup, the cluster will be kept in the list for future trials")
			outgroup_context.next_selection_must_have_outgroup = False
		else:
			# there are no candidat clsuters that might change the outgroup_to_concat_seq_ratio - remove it
			outgroup_context.add_cluster_to_not_selected(cid_to_check)
		return False
	if ratio_of_new_to_existing < new_to_existing_species_threshold_ratio and ratio_of_new_out_of_all < min_species_ratio_for_cluster:
		logger.debug(
			"Checking outgroup %s cluster_id %s - ratio_of_new_out_of_all is too small %f (smaller than %f and ratio_of_new_to_existing=%f) returning false" %
			(outgroup_context.outgroup_name, cid_to_check, ratio_of_new_out_of_all, min_species_ratio_for_cluster,
			 ratio_of_new_to_existing))
		outgroup_context.add_cluster_to_not_selected(cid_to_check)
		return False

	outgroup_context.add_cluster(cid_to_check)
	return True

def pick_next_cluster_for_common_NOoutgroup_greedy(ploidb_context, outgroup_context):
	logger.debug("Picking next cluster for %s" % outgroup_context.outgroup_name)
	cids_with_outgroup_seq = outgroup_context.candidate_cids_with_outgroup_seq
	cids_without_outgroup_seq = outgroup_context.candidate_cids_without_outgroup_seq
	logger.debug("cids_with_outgroup_seq=%s" % ",".join(cids_with_outgroup_seq))
	logger.debug("cids_without_outgroup_seq=%s" % ",".join(cids_without_outgroup_seq))

	cluster_wo_outgroup = None
	if len(cids_without_outgroup_seq) > 0 and not outgroup_context.next_selection_must_have_outgroup:
		next_cid = cids_without_outgroup_seq[0]
		cluster_wo_outgroup = ploidb_context.get_cluster_by_id(next_cid)
	logger.debug("next cluster without outgroup is %s " % (cluster_wo_outgroup))
	clusters_to_check = list()
	if cluster_wo_outgroup is not None: clusters_to_check.append(cluster_wo_outgroup)
	ordered_clusters = sorted(clusters_to_check, key=lambda cluster: cluster.get_data_matrix_size(estimated=True),
							  reverse=True)
	cluster_was_added = False
	for c in ordered_clusters:
		if cluster_was_added: break
		cluster_was_added = check_and_add_cluster_to_concat(c.cluster_id, ploidb_context, outgroup_context)


def pick_next_cluster_for_common_outgroup_greedy(ploidb_context, outgroup_context):
	logger.debug("Picking next cluster for %s" % outgroup_context.outgroup_name)
	cids_with_outgroup_seq = outgroup_context.candidate_cids_with_outgroup_seq
	cids_without_outgroup_seq = outgroup_context.candidate_cids_without_outgroup_seq
	logger.debug("cids_with_outgroup_seq=%s" % ",".join(cids_with_outgroup_seq))
	logger.debug("cids_without_outgroup_seq=%s" % ",".join(cids_without_outgroup_seq))
	cluster_w_outgroup = None
	#Check if there are cluster ids with Outgroup:
	if len(cids_with_outgroup_seq) > 0:
		next_cid = cids_with_outgroup_seq[0]
		cluster_w_outgroup = ploidb_context.get_cluster_by_id(next_cid)
	#Check if there are cluster ids without Outgroup:
	cluster_wo_outgroup = None
	if len(cids_without_outgroup_seq) > 0 and not outgroup_context.next_selection_must_have_outgroup:
		next_cid = cids_without_outgroup_seq[0]
		cluster_wo_outgroup = ploidb_context.get_cluster_by_id(next_cid)

	#Next candidates with and without outgroup:
	logger.debug("next cluster with outgroup is %s without %s" % (cluster_w_outgroup, cluster_wo_outgroup))
	clusters_to_check = list()
	if cluster_w_outgroup is not None: clusters_to_check.append(cluster_w_outgroup)
	if cluster_wo_outgroup is not None: clusters_to_check.append(cluster_wo_outgroup)
	ordered_clusters = sorted(clusters_to_check, key=lambda cluster: cluster.get_data_matrix_size(estimated=True),
							  reverse=True)
	cluster_was_added = False
	for c in ordered_clusters:
		if cluster_was_added: break
		cluster_was_added = check_and_add_cluster_to_concat(c.cluster_id, ploidb_context, outgroup_context)


def get_best_clusters_coverage_for_outgroup(conn, ploidb_context, outgroup_name, outgroup_type, num_of_clusters,
											data_matrix_size):
	logger.debug(
		"Processing the outgroup %s type=%s data_matrix_size=%d" % (outgroup_name, outgroup_type, data_matrix_size))

	clusters_ids_with_outgroup_seq, clusters_ids_without_outgroup_seq = get_clusters_lists_per_outgroup(conn,
																										ploidb_context,
																										outgroup_name)
	outgroup_context = OutGroupSelectionContext(ploidb_context, outgroup_name, outgroup_type,
												clusters_ids_with_outgroup_seq, clusters_ids_without_outgroup_seq)

	while outgroup_context.are_candidates_left():
		pick_next_cluster_for_common_outgroup_greedy(ploidb_context, outgroup_context)
		logger.debug("After next cluster pick %s" % outgroup_context)
	return outgroup_context


def get_best_clusters_coverage_for_outgroup_Mdebug(conn, ploidb_context, outgroup_name, outgroup_type, num_of_clusters,
											data_matrix_size):
	logger.debug(
		"Processing the outgroup %s type=%s data_matrix_size=%d" % (outgroup_name, outgroup_type, data_matrix_size))

	clusters_ids_with_outgroup_seq, clusters_ids_without_outgroup_seq = get_clusters_lists_per_outgroup(conn,
																										ploidb_context,
																										outgroup_name)
	outgroup_context = OutGroupSelectionContext(ploidb_context, outgroup_name, outgroup_type,
												clusters_ids_with_outgroup_seq, clusters_ids_without_outgroup_seq)
	while outgroup_context.are_candidates_left():
		pick_next_cluster_for_common_outgroup_greedy(ploidb_context, outgroup_context)
		logger.debug("After next cluster pick %s" % outgroup_context)
	return outgroup_context


def get_best_clusters_coverage_for_NOoutgroup(conn, ploidb_context, outgroup_name, outgroup_type, num_of_clusters,data_matrix_size):
	#logger.debug("Processing the outgroup %s type=%s data_matrix_size=%d" % (outgroup_name, outgroup_type, data_matrix_size))

	clusters_ids_with_outgroup_seq, clusters_ids_without_outgroup_seq = get_clusters_lists_per_outgroup(conn,ploidb_context,outgroup_name)

	outgroup_context = OutGroupSelectionContext(ploidb_context, outgroup_name, outgroup_type,clusters_ids_with_outgroup_seq, clusters_ids_without_outgroup_seq)
	logger.debug("clusters_ids_with_outgroup_seq")
	logger.debug(clusters_ids_with_outgroup_seq)
	logger.debug("clusters_ids_without_outgroup_seq")
	logger.debug(clusters_ids_without_outgroup_seq)

	while outgroup_context.are_candidates_left():
		pick_next_cluster_for_common_NOoutgroup_greedy(ploidb_context, outgroup_context)
		logger.debug("After next cluster pick %s" % outgroup_context)
	return outgroup_context


def Check_outgroup_forLegalName(outgroup_name):

	return 0 # MDebug: need to check what kind of check is needed for this version


def get_selected_outgroup(context_local):

	with open(context_local.outgroup_file) as f_outgroup:
		if context_local.UserFlags_dict['Outgroup_Flag'] == 'User':
			for line in f_outgroup:
				#outgroup_name = line.replace(' ','_').rstrip()
				outgroup_name = return_codedName(line)
				return outgroup_name
		else:
			for line in f_outgroup:
				if "(organism)" in line:
					outgroup_str = line[line.find("[")+1:line.find("]")]
					temp_out_name = outgroup_str.replace(" (organism)","")
					outgroup_name = temp_out_name.replace(" ","_").replace(':','_').replace("-","_").replace(",","").replace("'","").replace("(","_").replace(")","_") #in case of BOLD
					return outgroup_name
				elif "(genus)" in line:
					#print("Looking for outgroup...\n")
					outgroup_name = line[line.find("[")+1:line.find("]")]
					temp_out_1 = outgroup_name.split(" ")
					outgroup_name =temp_out_1[0]
					return outgroup_name
				elif "None" in line:
					return 'NULL'
	return 'FailedToRead_Outgroup'


def rerun_get_selected_outgroup(origin_outgroup_file):

	with open(origin_outgroup_file) as f_outgroup:
		for line in f_outgroup:
			if "(organism)" in line:
				outgroup_name = line[line.find("[")+1:line.find("]")]
				temp_out_1 = outgroup_name.split(" ")
				temp_out_2 = temp_out_1[0] + "_" + temp_out_1[1]
				temp_out_2 = temp_out_2.replace(':','_')    #in case of BOLD:3434 in name
				outgroup_name = temp_out_2
				return outgroup_name
			elif "(genus)" in line:
				#print("Looking for outgroup...\n")
				outgroup_name = line[line.find("[")+1:line.find("]")]
				temp_out_1 = outgroup_name.split(" ")
				outgroup_name =temp_out_1[0]
				return outgroup_name
			elif "None" in line:
				return 'NULL'
			else:
				#outgroup_name = line.replace(' ','_').rstrip()
				outgroup_name = return_codedName(line)
				return outgroup_name
	return 'FailedToRead_Outgroup'

def init_common_outgroup_MDebug(ploidb_context):
	# Open file with Outgroup outcome - either None or selected outgroup name:
	f_outgroup = open(ploidb_context.outgroup_file, 'w')
	# This way we ensure that the clusters table exists
	write_stats_to_db(ploidb_context)


	with get_db_conn(ploidb_context.db_file) as conn:
		curs = conn.cursor()

		# Iterating all possible outgroups - the WHERE condition is only optimization to make sure very small outgroups are not checked
		query = "select outgroup_name,outgroup_type,num_of_clusters,data_matrix_size " \
				"from outgroup_stats " \
				"where data_matrix_size >= 0.25 * (select max(data_matrix_size) from outgroup_stats) limit 250"

		curs.execute(query)
		outgroup_rows = curs.fetchall()
		logger.info("Found %d candidates for outgroup" % len(outgroup_rows))
		logger.info(outgroup_rows)
		#get all genera in current Alignment species list:
		genera_list=[]
		for specie in ploidb_context.species_list_names:
			specie_split=specie.split(' ')
			genera_list.append(specie_split[0])
			logger.debug("Init Common Outgroup, genera to exclude: %s\n" % specie_split[0])
		outgroup_context_list = list()
		get_next_outgroup=0
		for outgroup_row in outgroup_rows:
			outgroup_name = outgroup_row[0]
			#Check if legal outgroup
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
		#ploidb_context.cluster_contexts_for_concat_tree = ploidb_context.get_clusters_by_ids(
		#	selected_outgroup_context.get_selected_cids())
		ploidb_context.top_outgroup_contexts = top_outgroup_contexts
		ploidb_context.outgroupSelection = get_selected_outgroup(ploidb_context)
		with open(ploidb_context.summary_file,'a') as f_sum:
			f_sum.write("Selected Outgroup: %s\n" % ploidb_context.outgroupSelection)


		Write_Clusters_data_File(ploidb_context)
		write_outgroup_for_concat_to_file(ploidb_context, selected_outgroup_context)



def init_common_outgroup(ploidb_context):
	# Open file with Outgroup outcome - either None or selected outgroup name:
	f_outgroup = open(ploidb_context.outgroup_file, 'w')
	# This way we ensure that the clusters table exists
	write_stats_to_db(ploidb_context)


	with get_db_conn(ploidb_context.db_file) as conn:
		curs = conn.cursor()

		# Iterating all possible outgroups - the WHERE condition is only optimization to make sure very small outgroups are not checked
		query = "select outgroup_name,outgroup_type,num_of_clusters,data_matrix_size " \
				"from outgroup_stats " \
				"where data_matrix_size >= 0.25 * (select max(data_matrix_size) from outgroup_stats) limit 250"

		curs.execute(query)
		outgroup_rows = curs.fetchall()
		logger.info("Found %d candidates for outgroup" % len(outgroup_rows))
		logger.info(outgroup_rows)
		#get all genera in current Alignment species list:
		genera_list=[]
		for specie in ploidb_context.species_list_names:
			specie_split=specie.split(' ')
			genera_list.append(specie_split[0])
			logger.debug("Init Common Outgroup, genera to exclude: %s\n" % specie_split[0])
		outgroup_context_list = list()
		get_next_outgroup=0
		for outgroup_row in outgroup_rows:
			outgroup_name = outgroup_row[0]
			#Check if legal outgroup
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
		ploidb_context.outgroupSelection = get_selected_outgroup(ploidb_context)
		with open(ploidb_context.summary_file,'a') as f_sum:
			f_sum.write("Selected Outgroup: %s\n" % ploidb_context.outgroupSelection)

		for cluster_context in ploidb_context.cluster_contexts_for_concat_tree:
			cluster_context.is_used_for_concat_tree = True

		Write_Clusters_data_File(ploidb_context)
		write_outgroup_for_concat_to_file(ploidb_context, selected_outgroup_context)

#changed the query of outgroup to bring only names that don't contain genus names of alignment:


def init_NO_common_outgroup(ploidb_context,outgroup_selection):

	logger.info("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
	logger.info("^^^^^^^^^ outgroup_selection is %s ^^^^^^^^^^^^^^^" % outgroup_selection)
	logger.info("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")

	# This way we ensure that the clusters table exists
	write_stats_to_db(ploidb_context)
	ploidb_context.UserOutgroupInc = 'NO'

	with get_db_conn(ploidb_context.db_file) as conn:
		curs = conn.cursor()

		# Iterating all possible outgroups - the WHERE condition is only optimization to make sure very small outgroups are not checked
		query = "select outgroup_name,outgroup_type,num_of_clusters,data_matrix_size " \
				"from outgroup_stats " \
				"where data_matrix_size >= 0.25 * (select max(data_matrix_size) from outgroup_stats) limit 50"

		curs.execute(query)
		outgroup_rows = curs.fetchall()
		logger.info("Found %d candidates for outgroup" % len(outgroup_rows))
		#get all genera in current Alignment species list:
		genera_list=[]
		for specie in ploidb_context.species_list_names:
			specie_split=specie.split(' ')
			genera_list.append(specie_split[0])
			logger.debug("Init Common Outgroup, genera to exclude: %s\n" % specie_split[0])
		outgroup_context_list = list()
		get_next_outgroup=0

		#------------------------------------ Additional Code for Cluster Selection (michal) ---------------------------------------#
		#New Addition of user putgroup -> need to check if included in one of the selected clusters, otherwise notify the user:

		total_calc_length = 0
		total_species_num = len(ploidb_context.species_list_names)
		total_species_coverage_list = []
		max_concat_seq_length = int(ott_config['common outgroup']['max_concat_seq_length'])
		logger.debug("Maximum concat length is %s" % str(max_concat_seq_length))

		# Initialize list and params for concat clusters:
		Chosen_cluster_ids_list=[]
		max_cluster_index_Id_dict={}
		all_cluster_length = 0
		Chosen_cluster_list=[]
		ref_species_list=[]
		num_of_clusters = len(ploidb_context.cluster_contexts)
		i=0
		#This section will determine which clusters will be added to the final concat:
		#First it will add clusters that add the most data (according to wheight as defined below)
		#Then it will continue to add clusters with wheight 0 according to the NumOf taxa in them (max first)
		#it will stop when the sum of median length of clusters is less then 20000
		while all_cluster_length < max_concat_seq_length:# or i < num_of_clusters:
			added_species_dict={}
			weight_dict={}
			length_dict={}
			for cluster_context in ploidb_context.cluster_contexts:
				max_cluster_index_Id_dict[cluster_context.index]=cluster_context.cluster_id
				if cluster_context.index in Chosen_cluster_list:
					continue
				else:
					cluster_len = get_cluster_ids_total_seq_len(ploidb_context, [cluster_context.cluster_id])
					current_context = ploidb_context.get_cluster_by_id(cluster_context.cluster_id)
					logger.debug(cluster_context.index)
					#NEED to add filter for cluster selection in case of no outgroup:
					#def get_cluster_contribution(outgroup_context, ploidb_context, context_to_check):
					cluster_species_list = current_context.get_list_of_species()
					logger.debug("Checking if to add the following cluster %s, length = %s, species count = %s" % (ploidb_context.get_cluster_by_id(cluster_context.cluster_id),str(cluster_len),str(len(cluster_species_list))))
					logger.debug("Chosen Cluster List: %s" %Chosen_cluster_list)
					logger.debug("Species List: %s" %cluster_species_list)
					added_species = list(set(cluster_species_list) - set(ref_species_list))
					logger.debug("diff_species List: %s" %added_species)
					weight_dict[cluster_context.index]=len(cluster_species_list)*len(added_species)  #Weight is the number of species in cluster TIMES the number os New species it adds to the concat
					logger.debug("Weight %s" %str(weight_dict[cluster_context.index]))
					length_dict[cluster_context.index]=cluster_len
					added_species_dict[cluster_context.index] = added_species
					logger.debug(weight_dict)
					logger.debug(length_dict)
			i+=1
			max_cluster_index = max(weight_dict, key=weight_dict.get)
			all_cluster_length+=length_dict[max_cluster_index]
			Chosen_cluster_list.append(max_cluster_index)
			Chosen_cluster_ids_list.append(max_cluster_index_Id_dict[max_cluster_index])
			ref_species_list.extend(added_species_dict[max_cluster_index])
			logger.debug("Updated clusters list: %s" % Chosen_cluster_list)
			logger.debug("Added species list: %s" % added_species_dict[max_cluster_index])
			logger.debug("Adding Cluster %s, Total length: %s" %(str(max_cluster_index),str(all_cluster_length)))
			#Check user outgroup:
			if 'Outgroup_User' in ploidb_context.UserFlags_dict:
				for taxID_str in added_species_dict[max_cluster_index]:
					if taxID_str == ploidb_context.UserOutgroupTaxId:
						logger.debug("Found User Outgroup in added cluster !!!")
						ploidb_context.UserOutgroupInc = 'YES'
				#if str(ploidb_context.UserOutgroupTaxId) in added_species_dict[max_cluster_index]:
				#	logger.debug("Found User Outgroup in added cluster !!!")
				#	ploidb_context.UserOutgroupInc = 'YES'
			if len(Chosen_cluster_list) ==  num_of_clusters:
				break


		#Update summary file in case User outgroup is not included in the output:
		if 'Outgroup_User' in ploidb_context.UserFlags_dict and ploidb_context.UserOutgroupInc == 'NO':
			with open(ploidb_context.summary_file,'w') as sum_f:
				sum_f.write("User Outgroup (%s) is Not included in the final Alignment\n" %ploidb_context.UserFlags_dict['Outgroup_User'])
			ploidb_context.UserOutgroupInMSA = False
			#statusFail_LogFile(ploidb_context, 'The User Outgroup you specified was not included in the final MSA and so OneTwoTree did not reconstruct a phylogeny.')
			#raise Exception(
			#	"Failed - The User Outgroup you specified was not included in the final MSA and so OneTwoTree did not reconstruct a phylogeny.")
		elif 'Outgroup_User' in ploidb_context.UserFlags_dict:
			#Write user outgroup name to file
			with open(ploidb_context.outgroup_file, 'w') as out_f:
				out_f.write(ploidb_context.UserFlags_dict['Outgroup_User'])
			ploidb_context.outgroupSelection = get_selected_outgroup(ploidb_context)
			with open(ploidb_context.summary_file,'a') as f_sum:
				f_sum.write("User Selected Outgroup: %s\n" % ploidb_context.outgroupSelection)
		elif outgroup_selection == 'None':
			with open(ploidb_context.outgroup_file, 'w') as out_f:
				out_f.write('None')
		#Copy chosen clusters to final list:
		logger.debug("Chosen Cluster index: %s" % Chosen_cluster_list)
		logger.debug("Chosen Cluster Ids: %s " % Chosen_cluster_ids_list)
		ploidb_context.cluster_contexts_for_concat_tree = ploidb_context.get_clusters_by_ids(Chosen_cluster_ids_list)
		logger.debug(ploidb_context.cluster_contexts_for_concat_tree)

		#for cluster_context in ploidb_context.cluster_contexts:
		for cluster_context in ploidb_context.cluster_contexts_for_concat_tree:
			if cluster_context.index in Chosen_cluster_list:
				logger.debug("Adding Cluster %s to concat list" %str(cluster_context.index ))
				cluster_context.is_used_for_concat_tree = True


		#Update SummaryDir with all clusters data:
		Write_Clusters_data_File(ploidb_context)
		logger.info("Saving clusters data in %s" % ploidb_context.summary_clusters_data_file)

def init_NO_common_outgroup_MDebug(ploidb_context,outgroup_selection):

	logger.info("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
	logger.info("^^^^^^^^^ outgroup_selection is %s ^^^^^^^^^^^^^^^" % outgroup_selection)
	logger.info("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")

	# This way we ensure that the clusters table exists
	write_stats_to_db(ploidb_context)
	ploidb_context.UserOutgroupInc = 'NO'

	with get_db_conn(ploidb_context.db_file) as conn:
		curs = conn.cursor()

		# Iterating all possible outgroups - the WHERE condition is only optimization to make sure very small outgroups are not checked
		query = "select outgroup_name,outgroup_type,num_of_clusters,data_matrix_size " \
				"from outgroup_stats " \
				"where data_matrix_size >= 0.25 * (select max(data_matrix_size) from outgroup_stats) limit 50"

		curs.execute(query)
		outgroup_rows = curs.fetchall()
		logger.info("Found %d candidates for outgroup" % len(outgroup_rows))
		#get all genera in current Alignment species list:
		genera_list=[]
		for specie in ploidb_context.species_list_names:
			specie_split=specie.split(' ')
			genera_list.append(specie_split[0])
			logger.debug("Init Common Outgroup, genera to exclude: %s\n" % specie_split[0])
		outgroup_context_list = list()
		get_next_outgroup=0

		#Update SummaryDir with all clusters data:
		Write_Clusters_data_File(ploidb_context)
		logger.info("Saving clusters data in %s" % ploidb_context.summary_clusters_data_file)


def create_taxId_acc_list(cluster_fasta_seqs):

	taxId_Acc_list=[]
	with open(cluster_fasta_seqs,'r') as f_sequence:
		for line in f_sequence:
			if '>gi' in line:
				taxid_data_list = line.split('|')
				taxId = taxid_data_list[3]
				specie_name = taxid_data_list[5]
				logger.debug("****************** create_taxId_acc_list - name issue *************************")
				logger.debug(specie_name)
				name = specie_name.replace(",","").replace("'","").replace("-"," ").replace("("," ").replace(")"," ")   #Raxml doesn't like ( in the name :)
				logger.debug(name)
				#AccessionId = taxid_data_list[7]
				AccessionId = taxid_data_list[1]
				taxId_Acc_list.append(name + '|' + taxId + '|' + AccessionId)
	return taxId_Acc_list



def Write_Accession_Loci_Matrix(ploidb_context):

	f_accession_clusters = open(ploidb_context.summary_accession_matrix_file,'w')
	f_accession_clusters.write('"Specie\Loci",')
	clusters_header = ['ClusterID','Desc','Type','length','SpeciesCnt','SpeciesIds']

	Dict_TaxName=dict()
	ClusterId_list=[]
	TaxId_list=[]
	TaxId_Clust_Cluster_dictionary=dict()

	# 1. Need TaxID-Name dictionary and ClusterNum-Desc dictionary
	with open (ploidb_context.summary_clusters_data_file,'r') as clusters_f:
		logger.debug(clusters_f)
		dr = csv.DictReader(clusters_f, delimiter=',',fieldnames=['ClusterID','Desc','Type','length','SpeciesCnt','SpeciesIds'])
		for next_row in dr:
			if next_row['ClusterID'].isdigit():
				check_ = next_row['SpeciesIds'].replace('[','')
				check__ = check_.replace(']','').replace("', '",',')
				list_of_accession_data = (check__.strip("'")).split(",")
				print(list_of_accession_data)
				for item in list_of_accession_data:
					split_data=item.split('|') # 0 is name, 1 is tax, 2 is accession
					Name = split_data[0]
					TaxId = split_data[1]
					Acc = split_data[2]
					dble_key=()
					TaxId_Clust_Cluster_dictionary[(TaxId,next_row['ClusterID'])]=Acc
					Dict_TaxName[TaxId]=Name
					#Clust_Cluster_dictionary[next_row['ClusterID']].append(Acc)
					if next_row['ClusterID'] not in ClusterId_list:
						ClusterId_list.append(next_row['ClusterID'])
					if TaxId not in TaxId_list:
						TaxId_list.append(TaxId)
				f_accession_clusters.write('"%s-%s",' %(next_row['ClusterID'],next_row['Desc']))

		f_accession_clusters.write('\n')
		idx=0
		for TaxId in TaxId_list:
			f_accession_clusters.write('"%s-%s",' %(TaxId,Dict_TaxName[TaxId]))
			for clusterId in ClusterId_list:
				if (TaxId,clusterId) in TaxId_Clust_Cluster_dictionary.keys():
					print(TaxId_Clust_Cluster_dictionary[(TaxId,clusterId)])
					f_accession_clusters.write(TaxId_Clust_Cluster_dictionary[(TaxId,clusterId)])
				else:
					f_accession_clusters.write('"-"')
				f_accession_clusters.write(',')
			f_accession_clusters.write('\n')
	return


#---------------------------------------------------------------------------------------
# Create the summary file for all clusters selected for the concat file:
def Write_Clusters_data_File(ploidb_context):

	list_species_concat=[]
	f_summary_clusters = open(ploidb_context.summary_clusters_data_file,'w')
	f_summary_clusters.write("ClusterID,Desc,Type,length,SpeciesCnt,SpeciesIds,IncludedInConcat\n")
	cluster_json_dict={}
	for cluster_context in ploidb_context.cluster_contexts:
		current_context = ploidb_context.get_cluster_by_id(cluster_context.cluster_id)
		cluster_species_list = current_context.get_list_of_species()
		ClusterID = str(cluster_context.index)
		#Check if cluster was chosen for the concat file:
		logger.debug("ploidb_context.species_names_by_tax_id------------------------ REMOVE")
		logger.debug(ploidb_context.species_names_by_tax_id)
		if cluster_context in ploidb_context.cluster_contexts_for_concat_tree:
			IncludedInConcat = 'yes'
			#count species in concat:
			list_species_concat.extend(cluster_species_list)
			for taxID in cluster_species_list:
				if taxID not in ploidb_context.final_species_names_by_tax_id.keys():
					ploidb_context.final_species_names_by_tax_id[taxID] = ploidb_context.species_names_by_tax_id[taxID]
		else:
			IncludedInConcat = 'no'
		#if cluster_context.cluster_desc != 'None':
		if str(cluster_context.gene_name) != 'None':
			desc_str = cluster_context.gene_name
		else:
			desc_str = cluster_context.cluster_desc
		#create taxid-accession dictionary:
		taxid_accession_list = create_taxId_acc_list(cluster_context.all_seqs_fasta_filename_no_multiple_accessions)
		#desc_str=desc_str.split(' ', 2)[2]
		cluster_len = str(get_cluster_ids_total_seq_len(ploidb_context, [cluster_context.cluster_id]))
		cluster_data_set = [ClusterID , desc_str , cluster_len , len(cluster_species_list) , cluster_species_list , IncludedInConcat]
		cluster_json_dict[ClusterID] = cluster_data_set
		#Create dict to build accesion loci matrix file:
		ploidb_context.accession_loci_matrix_dict[ClusterID] = taxid_accession_list
		ploidb_context.loci_description_dict[ClusterID] = desc_str
		#f_summary_clusters.write(""""%s","%s","%s","%s","%d","%s","%s"\n""" %(ClusterID,desc_str,cluster_context.cluster_type_dict[ClusterID],cluster_len,len(cluster_species_list),cluster_species_list,IncludedInConcat))
		if IncludedInConcat == 'yes':
			f_summary_clusters.write(""""%s","%s","%s","%s","%d","%s","%s"\n""" %(ClusterID,desc_str,cluster_context.cluster_type_dict[ClusterID],cluster_len,len(cluster_species_list),taxid_accession_list,IncludedInConcat))
		else:
			#Clusters not included will be zipped together for download:
			logger.debug("Cluster %s not included in final alignment:" %ClusterID)
			logger.debug(""""%s","%s","%s","%s","%d","%s","%s"\n""" %(ClusterID,desc_str,cluster_context.cluster_type_dict[ClusterID],cluster_len,len(cluster_species_list),taxid_accession_list,IncludedInConcat))

	f_summary_clusters.close()
	#Write_Accession_Loci_Matrix(ploidb_context,ploidb_context.accession_loci_matrix_dict,ploidb_context.loci_description_dict,list(set(list_species_concat)))
	Write_Accession_Loci_Matrix(ploidb_context)

	with open(ploidb_context.summary_file,'a') as f_sum:
		#f_sum.write("Species with large data: %d\n" %len(ploidb_context.large_taxid_seqs_number_dict))
		f_sum.write("Species in final alignment file: %d\n" %len(set(list_species_concat)))
		#f_sum.write("TaxIds with no data in Genbank: %s\n" % len(ploidb_context.taxIds_not_in_genbank_list))
	with open(ploidb_context.summary_dir+'/Final_TaxId_vs_Name.txt','w') as f_final:
		for taxId in ploidb_context.final_species_names_by_tax_id.keys():
			f_final.write("%s,%s\n" %(taxId,ploidb_context.final_species_names_by_tax_id[taxId]))
	with open(ploidb_context.large_species_data,'w') as large_f:
		large_f.write("TaxId,NumberOfSeqs\n")
		for largeTax in ploidb_context.large_taxid_seqs_number_dict.keys():
			large_f.write("%s,%s\n" %(largeTax,ploidb_context.large_taxid_seqs_number_dict[largeTax]))
	#Create Constraint File if chosen by the user:
	#ploidb_context.UserFlags_dict["TaxaConstraint"] = 'None'#options are Genus / Family

	with open(ploidb_context.working_dir + '/ConstraintTree_user.txt') as const_file_h:
		if is_empty_file(const_file_h) == 0: #File is empty
			ploidb_context.UserFlags_dict["ConstraintTree"] = 'None'

	if ploidb_context.UserFlags_dict["ConstraintTree"] != 'None':
		# Set the var to notify this tree will be generated with a constraint file:
		ploidb_context.UserFlags_dict["Constraint_Flag"] = 'On'
		Constraint_Type = ploidb_context.UserFlags_dict["ConstraintTree"]
		#logger.debug("Create constraint tree according to user request: %s" % )
		if 'TaxaConstOp' in ploidb_context.UserFlags_dict.keys():
			Constraint_Taxa_type = ploidb_context.UserFlags_dict["TaxaConstOp"]
			logger.debug("Constraint taxa type: %s" % Constraint_Taxa_type)
		if 'UserConstOp' in ploidb_context.UserFlags_dict.keys():
			Constraint_User_type = ploidb_context.UserFlags_dict["UserConstOp"]
			logger.debug("Constraint user type: %s" % Constraint_User_type)
		#Call the relevant function to create the constraint tree file:
		if Constraint_Type == 'TaxaConst':	#Genus/Family/Both
			if Constraint_Taxa_type == 'GenusFamily_const':
				Create_Constraint_Both(ploidb_context.working_dir,ploidb_context.final_species_names_by_tax_id)
			else:
				create_ConstraintTree(ploidb_context.working_dir,ploidb_context.final_species_names_by_tax_id,Constraint_Taxa_type)
		elif  Constraint_Type == 'UserConst': #TaxID list/Constraint tree (need to remove species not included
			if Constraint_User_type == 'TaxaList':	#### NOT DONE !!!!
				create_ConstraintTree(ploidb_context.working_dir, ploidb_context.final_species_names_by_tax_id,
									  Constraint_User_type)
	return



def get_best_outgroup_context(outgroup_context_list):
	newlist = sorted(outgroup_context_list, key=lambda outgroup_context: outgroup_context.get_outgroup_score(),
					 reverse=True)
	max_top_outgroups = 40
	if len(newlist) < 40:
		max_top_outgroups = len(newlist)-1
	return newlist[0], newlist[0:max_top_outgroups]


def write_outgroup_for_concat_to_file(ploidb_context, outgroup_context):
	for cluster_context in outgroup_context.get_clusters_with_outgroup():
		out_seqs = get_common_outgroup_from_db_results(ploidb_context,ploidb_context.db_file, cluster_context.cluster_id,
													   None, None, outgroup_context.outgroup_name,
													   outgroup_context.outgroup_type)

		cluster_context.no_of_outgroup_sequence_for_concat = len(out_seqs)
		if len(out_seqs) == 0:
			logger.info(
				"While init common outgroup - no outgroup was found for cluster %s " % cluster_context.cluster_id)

		output_filename = cluster_context.outgroup_sequence_for_concat_filename
		out_handle = open(output_filename, 'w')
		logger.debug("Writing %i out group seqs the following to output file %s" % (len(out_seqs), output_filename))
		# todo: change the name so that it would be of the genus. maybe do it where all the names are being changed
		write_standard_seqs(out_handle, out_seqs)
		#SeqIO.write(out_seqs, out_handle, "fasta")
		out_handle.close()


def is_non_zero_file(fpath):
	return True if os.path.isfile(fpath) and os.path.getsize(fpath) > 0 else False
