import configparser
import hashlib
import io
import json
import os
import pickle
import re
import csv
from Bio import SeqIO, Phylo
from Bio import Entrez
import logging
import shutil
import subprocess
import sys
import sqlite3
import string
import unicodedata
from Bio.Align.Applications import MafftCommandline

# Because of long lines in clusters files and accession tables we enlarged the size to enable csv reader to read all lines:
csv.field_size_limit(2000000)

__author__ = 'moshe'

# TODO: these are already in the config file!!! kill this multiplication!!!!
#genbankDir = "/groups/itay_mayrose/mayaschuster/genebank/pln_2014_04"
#blastBinPath = "/groups/itay_mayrose/mosheein/ncbi-blast-2.2.26+/bin/"
#blastDbPath = "/biodb/BLAST/Nucleotides/nt"



#GROUP_LIST = ['rod','pri','mam','pln','vrl','vrt']

#GROUP_LIST = ['pri','rod','mam','vrt','inv','pln','vrl']
# #GROUP_LIST_FILES_dict = {'pri':56, 'rod':30, 'mam':39,'vrt':64,'inv':152,'pln':145,'vrl':48} #Updtae GenBank2017

GROUP_LIST = ['pri','rod','mam','vrt','inv','pln']
GROUP_LIST_FILES_dict = {'pri':57, 'rod':30, 'mam':39,'vrt':80,'inv':157,'pln':157}


#GROUP_LIST_FILES_dict = {'rod':30, 'pri':56, 'mam':39,'pln':145,'vrl':48,'vrt':64}


#GROUP_LIST = ['pln']  # Add ver
#GROUP_LIST_FILES_dict = {'pln':143}


#GROUP_LIST = ['rod','pri','mam','pln']  # Add ver
#GROUP_LIST_FILES_dict = {'rod':31, 'pri':53, 'mam':37,'pln':145}

ALIGN_UNUSED_CLUSTERS = True


logger = logging.getLogger('ploiDB-main')
ott_config = configparser.ConfigParser()



#running dos2unis on all files for concat
def run_dos2unix(index_fileList_f):
	#rename_indexed_files(index_fileList_f)
	logger.debug(" Running dos2unix on all msa's for concat file:")
	with open(index_fileList_f,'r') as indexed_list_f:
		for f_line in indexed_list_f:
			path_file = f_line.rstrip()
			logger.debug('dos2unix ' + path_file)
			os.system('dos2unix ' + path_file)
	return

def rename_indexed_files(index_fileList_f):
	logger.debug("Rename files to concat according to cluster:")
	shutil.copyfile(index_fileList_f,index_fileList_f+'_original')
	with open(index_fileList_f+'_original','r') as indexed_list_f:
		f_toWrite=open(index_fileList_f,'w')
		for line in indexed_list_f:
			r = re.compile('concat/(.*?)/')
			m = r.search(line)
			file_to_rename = line.rstrip()
			new_file_name = file_to_rename + "_" + m.group(1)
			os.rename(file_to_rename,new_file_name)
			f_toWrite.write(new_file_name+'\n')
	return


#
# Helper method for getting the number of a results that will be returned for a given Entrez query
#
def getResultsCountForDb_old(term, dbName):
	handle = Entrez.egquery(term=term)
	record = Entrez.read(handle)
	handle.close()

	for row in record["eGQueryResult"]:
		if row["DbName"] == dbName:
			print("Found " + row["Count"] + " results for db=" + dbName + " term=" + term)  # change this to logger.info
			return int(row["Count"])


def initCommandLineLogger():
	logger.setLevel(logging.DEBUG)
	# create console handler and set level to debug
	ch = logging.StreamHandler(sys.stdout)
	ch.setLevel(logging.DEBUG)
	# create formatter
	formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
	# add formatter to ch
	ch.setFormatter(formatter)
	# add ch to logger
	logger.addHandler(ch)


def init_ott_config(config_filename):
	global ott_config
	logger.info("Reading ploidb config from %s" % config_filename)
	ott_config.read(config_filename)
	logger.info("ott_config=%s" % ott_config)


#
# Run external process
#
def runExternalCommand(commandToExecute):
	logger.debug("Executing the following command - " + commandToExecute)
	p = subprocess.Popen(commandToExecute, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	for line in p.stdout.readlines():
		print(line)
	retval = p.wait()


#
# Run external process
#
def runExternalCommandOutToFile(commandToExecute, outFilename):
	logger.info("Executing the following command - " + commandToExecute)
	outFile = open(outFilename, "w")

	p = subprocess.Popen(commandToExecute, shell=True, stdout=outFile, stderr=subprocess.STDOUT)
	retval = p.wait()
	outFile.close()
	return retval


#
# Run external process
#
def exec_external_command_redirect_output(command_to_exec, outfile=None, errfile=None, cwd=None):
	out_handle = None
	err_handle = None
	if outfile is None:
		command_out = subprocess.PIPE
	else:
		out_handle = open(outfile, "w")
		command_out = out_handle

	if errfile is None:
		command_err = subprocess.PIPE
	else:
		err_handle = open(errfile, "w")
		command_err = err_handle

	logger.info("Executing the following command %s - out written to %s error written to %s", command_to_exec, outfile,
				errfile)
	logger.info("cwd=%s" % cwd)
	p = subprocess.Popen(command_to_exec, shell=True, cwd=cwd, stdout=command_out, stderr=command_err)
	# Added by Michal - Stuck in some cases (large input files for formatdb
	stdout, stderr = p.communicate()
	retval = p.wait()

	if retval == 0:
		logger.info("Execution return code was %i for command %s" % (retval, command_to_exec))
	else:
		logger.error("Execution return code was %i for command %s" % (retval, command_to_exec))
		raise Exception("Failed cmd")

	if out_handle is not None:
		out_handle.close()
	if err_handle is not None:
		err_handle.close()

	return retval


def getPropertyFromFastaSeqHeader(seqHeader, property):
	# print("getting prop from " + seqHeader)
	# gi is the first, so it needs special treatment (I didn't have the strength writing a smart regexp)
	if property == "gi":
		exp = r"gi\|(.*?)(\||\\n)"
		groupIdx = 1
	else:
		exp = r"(\>|\|)%s\|(.*?)(\||\\n|$)" % property
		groupIdx = 2
	m = re.search(exp, seqHeader)

	if m is None:
		return None
	else:
		propertyValue = m.group(groupIdx)
		return propertyValue


def getTaxid(seq):
	speciesTaxId = None
	for feature in seq.features:
		if feature.type == "source":
			# TODO: db_xref doesn't always contain the taxonid
			db_xrefs = feature.qualifiers.get("db_xref")
			for db_xref in db_xrefs:
				parts = db_xref.split(":")
				if parts[0].strip() == "taxon":
					# print "looking for " + db_xref + " found " + parts[1]
					speciesTaxId = parts[1].strip()
					break
				else:
					speciesTaxId = None
	return speciesTaxId




''' no need for it - return everything in get_seq_organism
def get_seq_taxonomy_last(seq,synonym_dict):
	taxonomy_last = seq.annotations["taxonomy"][-1]
	organism = seq.annotations["organism"]
	logger.debug("Taxonomy Last: %s Organism: %s ",taxonomy_last, organism)
	# If the species of the sequence have a synonym name that's in the Genus - return it instead of the genebank name
	for syn in synonym_dict.keys():
		if organism in syn:
			logger.debug("Found synonym genus name for %s, changing it %s",taxonomy_last,synonym_dict[syn].split(" ")[0])
			taxonomy_last = synonym_dict[syn].split(" ")[0]
	return taxonomy_last
'''


def get_seq_organism(seq, synonym_dict):
	organism = seq.annotations["organism"]
	# logger.debug("Organism: %s",organism)
	# If the species of the sequence have a synonym name that's in the Genus - return it instead of the genebank name
	for syn in synonym_dict.keys():
		if organism == syn:
			logger.debug("Found synonym organism name for %s, changing it to %s", organism, synonym_dict[syn])
			organism = synonym_dict[syn]
	return organism


def restandarize_seq_short(seq):
	#seqid = getPropertyFromFastaSeqHeader(seq.description, "seqid")
	seqid = getPropertyFromFastaSeqHeader(seq.description, "gi")
	seq.description = "%s" % (seqid)


def restandarize_seq(seq, gi=None, taxonid=None, organism=None, original_organism=None, description=None, seqid=None,
					 original_taxonid=None, original_subsp_organism=None, original_subsp_taxonid=None):
	if gi is None:
		gi = getPropertyFromFastaSeqHeader(seq.description, "gi")
	if taxonid is None:
		taxonid = getPropertyFromFastaSeqHeader(seq.description, "taxonid")
	if organism is None:
		organism = getPropertyFromFastaSeqHeader(seq.description, "organism")

	seq_description = "gi|%s|taxonid|%s|organism|%s" % (
		gi,
		taxonid,
		organism)

	if original_taxonid is None:
		original_taxonid = getPropertyFromFastaSeqHeader(seq.description, "original_taxonid")
	if original_taxonid is not None:
		seq_description = "%s|original_taxonid|%s" % (
			seq_description,
			original_taxonid)

	if original_organism is None:
		original_organism = getPropertyFromFastaSeqHeader(seq.description, "original_organism")
	if original_organism is not None:
		seq_description = "%s|original_organism|%s" % (
			seq_description,
			original_organism)

	if original_subsp_taxonid is None:
		original_subsp_taxonid = getPropertyFromFastaSeqHeader(seq.description, "original_subsp_taxonid")
	if original_subsp_taxonid is not None:
		seq_description = "%s|original_subsp_taxonid|%s" % (
			seq_description,
			original_subsp_taxonid)

	if original_subsp_organism is None:
		original_subsp_organism = getPropertyFromFastaSeqHeader(seq.description, "original_subsp_organism")
	if original_subsp_organism is not None:
		seq_description = "%s|original_subsp_organism|%s" % (
			seq_description,
			original_subsp_organism)

	if description is None:
		description = getPropertyFromFastaSeqHeader(seq.description, "description")
	if seqid is None:
		#seqid = getPropertyFromFastaSeqHeader(seq.description, "seqid")
		seqid = getPropertyFromFastaSeqHeader(seq.description, "gi")

	seq.description = "%s|seqid|%s|description|%s" % (
		seq_description,
		seqid,
		description)


def standarizeSequence(seq, genus_name=None, override_organism=None, override_taxonid=None, original_organism=None):
	if genus_name is None:
		# TODO: use "ORGANISM" in the DB format to get the genus
		final_genus_name = "Unknown Genus"
	else:
		final_genus_name = genus_name

	# for ann in seq.annotations:
	#    print(ann)
	if override_organism is not None:
		organism = override_organism
	else:
		organism = seq.annotations["organism"]

	# taxonomy = seq.annotations["taxonomy"]
	#gi = seq.annotations["gi"]

	if override_taxonid is not None:
		taxon_id = override_taxonid
	else:
		taxon_id = getTaxid(seq)

	sequence_version = seq.annotations["sequence_version"]
	#print("gi=" + gi)
	#print("organism=" + organism)
	#print("sequence_version=" + str(sequence_version))
	#print("accessions[0]=" +  seq.annotations["accessions"][0])

	if original_organism is None:
		original_organism = organism

	seq.description = "gi|%s|taxonid|%s|organism|%s|genus|%s|original_organism|%s|seqid|%s|description|%s" % (
		#gi,
		seq.id,
		taxon_id,
		organism,
		final_genus_name,
		original_organism,
		seq.id,
		seq.description)


def write_standard_seq(out, seq):
	out.write(">%s\n%s\n" % (seq.description, seq.seq))


def write_standard_seqs(out, seqs):
	for seq in seqs:
		out.write(">%s\n%s\n" % (seq.description, seq.seq))


def getSpeciesSequenceFileLocation(taxonId):
	seqPerSpeciesDir = ott_config['general']['seq_per_species_dir']
	path = seqPerSpeciesDir + taxonId[0] + "/" + taxonId
	if os.path.exists(path):
		return path
	else:
		return None

def checkIfSpeciesExists(taxonId):
	genbank_db = ott_config['general']['DB_DIR'] + ott_config['concat_headir']['GENBANK_DB']
	grp_list = ott_config['general']['GENBANK_GRP_LIST'].split(',')
	conn = sqlite3.connect(genbank_db)
	out_query = ''
	curs = conn.cursor()
	count_grp = 1
	#Check relevant groups and search all tables accordingly:
	logger.debug("In checkIfSpeciesExists: check for taxID: %s" %taxonId)
	for grp_name in grp_list:
		grp_db_table = 'Genbank_seqs_' + grp_name
		out_query += "SELECT TaxonId from %s where TaxonID like '%s'" % (grp_db_table,taxonId)
		if count_grp < len(grp_list):
			out_query += " UNION "
			count_grp+=1
	logger.debug("Query is: %s" % out_query)
	curs.execute(out_query)
	rows = curs.fetchall()
	if taxonId in rows:
		conn.close()
		logger.debug("Found TaxId")
		return taxonId
	else:
		conn.close()
		logger.debug("No TaxId")
		return 'None'


def writeSequenceInFastaFormat(out, nextTaxId, seq):
	Debug_Flag = 0
	# for ann in seq.annotations:
	#    print(ann)
	organism = seq.annotations["organism"]
	# taxonomy = seq.annotations["taxonomy"]

	#gi = seq.annotations["gi"]
	if Debug_Flag==1: logger.debug(seq)
	if Debug_Flag==1: logger.debug(seq.annotations.keys())
	gi = seq.annotations["accessions"][0] #accessions

	sequence_version = seq.annotations["sequence_version"]
	out.write(">gi|%s|taxonid|%s|organism|%s|seqid|%s|description|%s\n%s\n" % (
		#gi,
		seq.id,
		nextTaxId,
		organism,
		seq.id,
		seq.description,
		seq.seq))


def get_ordered_seqs(seq_list):
	ordered_seqs = sorted(seq_list, key=lambda seq: seq.description)
	return ordered_seqs


def write_standarize_fasta(out_fasta_filename, seq_list):
	ordered_seqs = get_ordered_seqs(seq_list)
	with open(out_fasta_filename, 'w') as out_handle:
		for seq in ordered_seqs:
			seq.id=seq.description
			SeqIO.write(seq, out_handle, "fasta")
			#out_handle.write(">%s\n%s\n" % (seq.description, seq.seq))


def write_fasta_short_desc(out_fasta_filename, seq_list):
	with open(out_fasta_filename, 'w') as out_handle:
		for seq in seq_list:
			#new_seq_desc = getPropertyFromFastaSeqHeader(seq.description, "seqid")
			new_seq_desc = getPropertyFromFastaSeqHeader(seq.description, "gi")
			out_handle.write(">%s\n%s\n" % (new_seq_desc, seq.seq))


def writeSequenceInFastaFormatNoTaxid(out, genusName, seq):
	writeSequenceInFastaFormat(out, genusName, getTaxid(seq), seq)


def concat_fasta_cluster(clusters, workdir, concat_fasta, concat_outgroup_fasta=None):
	logger.info("Starting concatenation in of %d clusters to %s. outgroup will be written to %s" % (
		len(clusters), concat_fasta, concat_outgroup_fasta))
	create_dir_if_not_exists(workdir)
	indexFileName = workdir + "/index.txt"

	all_outgroup_recs = list()
	with open(indexFileName, "w") as handle:
		for cluster_context in clusters:
			aligned_fasta_filname = cluster_context.all_seqs_for_concat_tree_creation_fasta_filename
			logger.debug("Writing %s to %s" % (aligned_fasta_filname, indexFileName))
			handle.write(aligned_fasta_filname + "\n")
			if cluster_context.no_of_outgroup_sequence_for_concat > 0:
				outgroup_recs = list(
					SeqIO.parse(cluster_context.outgroup_sequence_for_concat_short_desc_filename, "fasta"))
				all_outgroup_recs.extend(outgroup_recs)

	unique_outgroup_recs = list()
	outgroup_names_recs = list()
	for outgroup_rec in all_outgroup_recs:
		if not outgroup_rec.name in outgroup_names_recs:
			logger.debug("Found unique outgroup %s" % outgroup_rec.name)
			unique_outgroup_recs.append(outgroup_rec)
			outgroup_names_recs.append(outgroup_rec.name)

	logger.debug("Found a total of %d otugroup seqs" % len(all_outgroup_recs))
	if concat_outgroup_fasta is not None:
		with open(concat_outgroup_fasta, "w") as concat_outgroup_handle:
			SeqIO.write(unique_outgroup_recs, concat_outgroup_handle, "fasta")

	#run_dos2unix(indexFileName)
	concat_align_command = "perl /groups/pupko/haim/pupkoSVN/trunk/programs/indelReliability/ConcateAlignments.pl %s %s %s NO NA NA" % \
						   (indexFileName, concat_fasta, workdir + "/concat-report.txt")
	exec_external_command_redirect_output(concat_align_command, workdir + "/concat.out", workdir + "/concat.err")


def createSeqCodesForFastaFiles(inferredorganismNameToOriginalSeqDescriptionList, outSeqCodesFilename):
	outSeqCodesHandle = open(outSeqCodesFilename, 'w')

	for inferredorganismNameToOriginalSeqDescription in inferredorganismNameToOriginalSeqDescriptionList:
		for organismName, OriginalSeqDescription in inferredorganismNameToOriginalSeqDescription.items():
			# outSeqCodesHandle.write(organismName + "\t" + OriginalSeqDescription + "\n")
			outSeqCodesHandle.write(OriginalSeqDescription + "\t" + organismName + "\n")
	outSeqCodesHandle.close()


def get_replaced_extension_filename(filename, old_ext, new_ext):
	return filename.replace(old_ext, new_ext)


def escape_organism_name(organism_name):
	escaped_organism = organism_name

	escaped_organism = escaped_organism.replace(":", "_") #BOLD:345345 case
	escaped_organism = escaped_organism.replace(" ", "_")
	escaped_organism = escaped_organism.replace(",", "_")
	escaped_organism = escaped_organism.replace("-", "_")
	escaped_organism = escaped_organism.replace("(", "")
	escaped_organism = escaped_organism.replace(")", "")

	escaped_organism = escaped_organism.replace("'", "_")
	escaped_organism = escaped_organism.replace("/", "_")
	escaped_organism = escaped_organism.replace("&", "AND")
	escaped_organism = unicodedata.normalize('NFKD', escaped_organism).encode('ascii', 'ignore')
	escaped_organism = str(escaped_organism, encoding='ascii')
	return escaped_organism


def rewrite_fasta_with_organism_name_as_desc(in_fasta_filename, out_organism_desc_fasta_filename, add_gi=False):
	out_handle = open(out_organism_desc_fasta_filename, 'w')

	names_dict = dict()
	for seq_record in SeqIO.parse(in_fasta_filename, "fasta"):
		organism = getPropertyFromFastaSeqHeader(seq_record.description, "organism")
		original_organism = organism

		logger.debug("rewriting %s" % organism)
		organism = escape_organism_name(original_organism)
		if add_gi:
			logger.debug("Adding GI=true was called")
			gi = getPropertyFromFastaSeqHeader(seq_record.description, "gi")
			out_handle.write(">%s_%s\n%s\n" % (organism, gi, seq_record.seq))
			logger.debug("rewrote to %s_%s" % (organism, gi))
		else:
			out_handle.write(">%s\n%s\n" % (organism, seq_record.seq))
			#To prevent names with Parentheses to appear in the final alignment, my cause Raxml to fail
			#organism.replace("(","").replace(")","")  added to escape_organism_name
			logger.debug("rewrote to %s" % organism)

		names_dict[original_organism] = organism
	return names_dict


def template(filename, vars):
	filein = open(filename)
	src = string.Template(filein.read())
	return (src.substitute(vars))


def rm_and_create_dir(dir):
	logger.debug("deleting and recreating dir %s" % dir)
	if os.path.exists(dir):
		shutil.rmtree(dir)
	os.makedirs(dir)


def create_dir_if_not_exists(dir):
	if not os.path.exists(dir):
		logger.debug("Creating %s" % dir)
		os.makedirs(dir)


def get_indexed_seq_hash(genbank_files_dir):
	files = []
	for grp_name in GROUP_LIST:
		number_of_genbank_files = GROUP_LIST_FILES_dict[grp_name]
		#files = [genbank_files_dir + "/gbpln%i.seq" % (i + 1) for i in range(number_of_genbank_files)]
		#genbank_files_dir = ott_config['general'][genbankDir + grp_name + '/']
		if not files:
			files = [genbank_files_dir + grp_name + "/gb"+grp_name+"%i.seq" % (i + 1) for i in range(number_of_genbank_files)]
		else:
			files.extend([genbank_files_dir + grp_name + "/gb"+grp_name+"%i.seq" % (i + 1) for i in range(number_of_genbank_files)])

	indexFilename = ott_config['general']['index_Filename']
	gb_vrl = SeqIO.index_db(indexFilename, files, "genbank")
	logger.debug("%i sequences indexed" % len(gb_vrl))
	return gb_vrl


def return_codedName(origin_name):

	tempName= origin_name.rstrip()
	tempName=tempName.replace(' ','_')
	CodedName = tempName.replace(':','_') #in case of a BOLD name
	return CodedName

def get_db_conn(db_file):
	conn = sqlite3.connect(db_file)
	return conn


def alter_filename(path, name_suffix):
	head, tail = os.path.split(path)
	(filename, ext) = os.path.splitext(tail)
	new_filename = os.path.join(head, filename + name_suffix + ext)
	return new_filename


def get_newick_tree_str(filename):
	tree = Phylo.parse(filename, 'newick').next()
	# output = StringIO.StringIO()
	output = io.StringIO()
	Phylo.write(tree, output, 'newick')
	return output.getvalue()


def get_out_filename_from_mb_config(mb_config_filename):
	mcmc_line_pattern = r"(\s*mcmc.*file=)(.*)\s*.*;"
	with open(mb_config_filename, "r") as mb_conf_handle:
		for line in mb_conf_handle:
			m = re.match(mcmc_line_pattern, line)
			if m:
				out_file = m.group(2)
				return out_file

	return None


def get_value_by_percentile(all_vals, percentile):
	all_vals_sorted = list(all_vals)
	all_vals_sorted.sort()

	percentiles = list()
	for percent in range(0, 100, 1):
		index_of_percentile = int(len(all_vals_sorted) * percent / 100)
		# logger.debug("Adding to percentiles index %d the value %d (which is in index %d in the original list)" % (percent,all_vals_sorted[index_of_percentile],index_of_percentile))
		percentiles.append(all_vals_sorted[index_of_percentile])

	index_of_percentile_val = int(len(all_vals_sorted) * percentile)

	return all_vals_sorted[index_of_percentile_val]


def write_stats_to_db(context):
	logger.info("Writing stats to DB")
	conn = sqlite3.connect(context.db_file)
	try:
		curs = conn.cursor()

		# Inserting genus data
		values_for_db = [(context.genus_id, context.id, len(context.cluster_contexts),
						  len(context.cluster_contexts_for_concat_tree),
						  len(context.cluster_contexts_for_loci_trees), len(context.unresolved_species_names))]
		logger.debug(values_for_db)
		query = "INSERT or REPLACE INTO genus_stats " \
				"(id, genus_name,num_clusters_total,num_clusters_concat,num_clusters_locus_tree,num_filtered_by_name) " \
				"VALUES (?, ?,?,?,?,?);"
		logger.debug("About to execute the following query %s" % query)
		curs.executemany(query, values_for_db)

		# Inserting cluster's data
		query = "INSERT or REPLACE INTO clusters " \
				"(id, genus_id,cluster_index,cluster_type,is_used_for_locus_tree,is_used_for_concat_tree,is_its,is_chloroplast,chloroplast_loci,num_of_taxons,align_length,data_matrix_size) " \
				"VALUES (?,?,?,?,?,?,?,?,?,?,?,?);"
		for cluster_context in context.cluster_contexts:
			c = cluster_context
			# c.init_data_matrix_size()
			values_for_db = [
				(cluster_context.cluster_id, context.genus_id, cluster_context.index, cluster_context.cluster_type,
				 cluster_context.is_used_for_locus_tree, cluster_context.is_used_for_concat_tree,
				 cluster_context.is_its_cluster, 0, 0, c.get_number_of_different_species(), c.get_cluster_length(True),
				 c.get_data_matrix_size(True))]
			logger.debug("About to execute the following query %s" % query)
			logger.debug(values_for_db)
			curs.executemany(query, values_for_db)

		conn.commit()
		# Inserting chloroplast data
		if context.chloroplast_cluster is not None:
			chloroplast_context = context.chloroplast_cluster
			c = chloroplast_context
			values_for_db = [(chloroplast_context.cluster_id, context.genus_id, chloroplast_context.index,
							  chloroplast_context.cluster_type,
							  chloroplast_context.is_used_for_locus_tree, chloroplast_context.is_used_for_concat_tree,
							  chloroplast_context.is_its_cluster, 1, len(chloroplast_context.cluster_list),
							  c.get_number_of_different_species(), c.get_cluster_length(True),
							  c.get_data_matrix_size(True))]
			logger.debug(values_for_db)
			logger.debug("About to execute the following query %s" % query)
			curs.executemany(query, values_for_db)

		conn.commit()
		# Inserting Mt data
		if context.Mt_cluster is not None:
			Mt_context = context.Mt_cluster
			c = Mt_context
			values_for_db = [(Mt_context.cluster_id, context.genus_id, Mt_context.index,
							  Mt_context.cluster_type,
							  Mt_context.is_used_for_locus_tree, Mt_context.is_used_for_concat_tree,
							  Mt_context.is_its_cluster, 1, len(Mt_context.cluster_list),
							  c.get_number_of_different_species(), c.get_cluster_length(True),
							  c.get_data_matrix_size(True))]
			logger.debug(values_for_db)
			logger.debug("About to execute the following query %s" % query)
			curs.executemany(query, values_for_db)

		conn.commit()
		# Inserting Nuc data
		if context.Nuc_cluster is not None:
			Nuc_context = context.Nuc_cluster
			c = Nuc_context
			values_for_db = [(Nuc_context.cluster_id, context.genus_id, Nuc_context.index,
							  Nuc_context.cluster_type,
							  Nuc_context.is_used_for_locus_tree, Nuc_context.is_used_for_concat_tree,
							  Nuc_context.is_its_cluster, 1, len(Nuc_context.cluster_list),
							  c.get_number_of_different_species(), c.get_cluster_length(True),
							  c.get_data_matrix_size(True))]
			logger.debug(values_for_db)
			logger.debug("About to execute the following query %s" % query)
			curs.executemany(query, values_for_db)

		conn.commit()
	finally:
		if conn is not None:
			conn.close()

	logger.info("DONE writing stats to DB")


def get_cluster_ids_total_data_matrix(ploidb_context, cluster_ids):
	total_data_matrix = 0
	for c in ploidb_context.get_clusters_by_ids(cluster_ids):
		total_data_matrix += c.get_data_matrix_size(True)
	return total_data_matrix


def hashfile(file_list, hasher, blocksize=65536):
	if not isinstance(file_list, list):
		logger.debug("Hashing based on %s" % file_list)
		file_list = [file_list]
	for file in file_list:
		logger.debug("Hashing based on %s" % ",".join(file_list))
		with  open(file, 'rb') as file_handle:
			buf = file_handle.read(blocksize)
			while len(buf) > 0:
				hasher.update(buf)
				buf = file_handle.read(blocksize)
	return hasher.hexdigest()


def get_number_of_different_species(fasta_file):
	return len(get_list_of_species(fasta_file))


def get_list_of_species(fasta_file):
	different_taxons = set()
	all_seqs = list(SeqIO.parse(fasta_file, "fasta"))
	for seq in all_seqs:
		taxonid = getPropertyFromFastaSeqHeader(seq.description, "taxonid")
		if taxonid not in different_taxons:
			different_taxons.add(taxonid)
	# logger.debug("Found %d different taxons for %d seqs in the file %s" % (len(different_taxons),len(all_seqs),fasta_file))
	return different_taxons


def get_cluster_ids_total_seq_len(ploidb_context, cluster_ids):
	# logger.debug("Getting total seq len for %s" % ",".join(cluster_ids))
	total_seq_len = 0
	for c in ploidb_context.get_clusters_by_ids(cluster_ids):
		#logger.debug("Seq length of %s is %d" % (c.cluster_id,c.seq_length))
		total_seq_len += c.get_cluster_length(estimated=True)
	return total_seq_len


def use_cache_results_if_exists(file_for_hash, file_to_check_in_cache, genus_name, cache_prefix):
	cache_flag = ott_config.getboolean('general', 'use_cache')

	if not cache_flag:
		logger.debug("use_cache flag is down - will not use cache ")
		return False

	if genus_name is None:
		logger.debug("genus_name is None - will not use cache ")
		return False

	hash = hashfile(file_for_hash, hashlib.md5())
	logger.debug("Hash returned is %s" % hash)

	cached_filename = os.path.join(ott_config['general']['cache_path'], genus_name)
	cached_filename = os.path.join(cached_filename, cache_prefix + hash)
	logger.debug("cached_filename = %s" % cached_filename)

	if not os.path.exists(cached_filename):
		logger.debug("No cache was found for %s using this filename %s" % (file_for_hash, cached_filename))
		return False

	if isinstance(file_to_check_in_cache, list):
		logger.debug("Using list of cached file %s" % ",".join(file_to_check_in_cache))
		# Cleaning the old stuff (not from cache)
		for next_file_to_check_in_cache in file_to_check_in_cache:
			logger.debug("Handling %s" % next_file_to_check_in_cache)
			if os.path.exists(next_file_to_check_in_cache) and os.path.isdir(next_file_to_check_in_cache):
				logger.debug("Deleting %s" % next_file_to_check_in_cache)
				shutil.rmtree(next_file_to_check_in_cache)
		# Verifying that cache of list, all the cache files/dirs should be on the same base dir
		target_paths = set()
		for next_file_to_check_in_cache in file_to_check_in_cache:
			base_path, dirname = os.path.split(next_file_to_check_in_cache)
			target_paths.add(base_path)
		if len(target_paths) != 1:
			raise Exception("Caching is only possible for the same target directory. Got %d directories for %s" % (
				len(target_paths), ",".join(file_to_check_in_cache)))
		else:
			final_base_path = list(target_paths)[0]
			logger.debug("All cahced files will be copied to %s" % final_base_path)
			my_copytree(cached_filename, final_base_path)
	else:
		if os.path.isdir(cached_filename):
			if os.path.exists(file_to_check_in_cache):
				shutil.rmtree(file_to_check_in_cache)
			shutil.copytree(cached_filename, file_to_check_in_cache)
		else:
			shutil.copy(cached_filename, file_to_check_in_cache)
		logger.info("Using cached version of %s for %s" % (cached_filename, file_to_check_in_cache))
	logger.debug("use_cache_results_if_exists - Returning true")
	return True


def my_copytree(src, dst, symlinks=False, ignore=None):
	logger.debug("Copying from %s to %s" % (src, dst))
	for item in os.listdir(src):
		s = os.path.join(src, item)
		d = os.path.join(dst, item)
		logger.debug("Copying %s to %s" % (s, d))
		if os.path.isdir(s):
			logger.debug("Copying using copytree")
			shutil.copytree(s, d, symlinks, ignore)
		else:
			logger.debug("Copying using copy2")
			shutil.copy2(s, d)


def copy_dir_to_another_dir(source_dir, target_dir):
	path, dirname = os.path.split(source_dir)
	target_dirname = os.path.join(target_dir, dirname)
	logger.debug("About to copy dir %s to %s", source_dir, target_dirname)
	shutil.copytree(source_dir, target_dirname)


def cache_file(file_for_hash, source_to_cache, genus_name, cache_prefix):
	cache_flag = ott_config.getboolean('general', 'use_cache')

	if not cache_flag:
		logger.debug("use_cache flag is down - will not use cache ")
		return False

	if genus_name is None:
		logger.debug("genus_name is None - will not cache any file")
		return

	hash = hashfile(file_for_hash, hashlib.md5())
	logger.debug("Hash returned is %s" % hash)

	cached_dir = os.path.join(ott_config['general']['cache_path'], genus_name)
	if not os.path.exists(cached_dir):
		logger.debug("Creating directory %s" % cached_dir)
		os.makedirs(cached_dir)

	cached_filename = os.path.join(cached_dir, cache_prefix + hash)
	logger.debug("Putting the file %s in cache as %s" % (source_to_cache, cached_filename))
	if os.path.exists(cached_filename) and os.path.isdir(cached_filename):
		logger.debug("Deleting existing cache at %s" % cached_filename)
		shutil.rmtree(cached_filename)

	if isinstance(source_to_cache, list):
		create_dir_if_not_exists(cached_filename)
		for next_file_or_folder_to_cache in source_to_cache:
			if os.path.isdir(next_file_or_folder_to_cache):
				copy_dir_to_another_dir(next_file_or_folder_to_cache, cached_filename)
			else:
				if not (os.path.exists(next_file_or_folder_to_cache)):
					logger.debug("Path %s NOT FOUND!!! Can't cache to ", next_file_or_folder_to_cache, cached_filename)
				else:
					logger.debug("About to copy %s to %s", next_file_or_folder_to_cache, cached_filename)
					shutil.copy(next_file_or_folder_to_cache, cached_filename)
	else:
		if os.path.isdir(source_to_cache):
			copy_dir_to_another_dir(source_to_cache, cached_filename)
		else:
			logger.debug("About to copy %s to %s", source_to_cache, cached_filename)
			shutil.copy(source_to_cache, cached_filename)

	# TODO: remove. this is for debugging the hash mechanism
	hash_filename = os.path.join(cached_dir, "src_" + cache_prefix + hash)
	shutil.copy(file_for_hash, hash_filename)


def exec_blast_all_vs_all(all_seqs_fasta_filename, fasta_file_for_blast, blast_results_filename, genus_name=None):
	logger.debug("About to execute blast all-vs-all on %s. new fasta will be %s" % (
		all_seqs_fasta_filename, fasta_file_for_blast))
	new_seqs = list()
	all_seqs = list(SeqIO.parse(all_seqs_fasta_filename, "fasta"))
	for seq in all_seqs:
		restandarize_seq_short(seq)
		new_seqs.append(seq)
	write_standarize_fasta(fasta_file_for_blast, new_seqs)

	# If we already have the blast all-vs-all results from previous executions - use it. Else run it
	get_cuu_dir = os.path.dirname(fasta_file_for_blast)
	logger.debug("Check if this is ITS cluster: %s" %fasta_file_for_blast)
	if not use_cache_results_if_exists(fasta_file_for_blast, blast_results_filename, genus_name, "ava"):
		logger.info("Preparing to running BLAST all-vs-all - calling formatdb")
		#fomat_db_command = "formatdb -i %s -pF" % fasta_file_for_blast
		fomat_db_command = "formatdb -i %s -pF" % fasta_file_for_blast
		exec_external_command_redirect_output(fomat_db_command)

		logger.info("Running BLAST all-vs-all")
		blast_command = "blastall -p blastn -d %s -i %s -v 100000 -b 100000 -e 1e-5 -m 8" % (
			fasta_file_for_blast, fasta_file_for_blast)
		exec_external_command_redirect_output(blast_command, blast_results_filename)
		cache_file(fasta_file_for_blast, blast_results_filename, genus_name, "ava")


def escape_str_for_sql_query(str):
	return str.replace("'", "''")


def get_first_seq_from_file(fasta_filename):
	seq = None
	if os.path.exists(fasta_filename):
		combined_recs = list(SeqIO.parse(fasta_filename, "fasta"))
		if len(combined_recs) > 0:
			seq = combined_recs[0]
	return seq


def restandarize_description_after_mafft(fasta_filename):
	seqs = list(SeqIO.parse(fasta_filename, "fasta"))
	for seq in seqs:
		if "_R_gi" in seq.description:
			seq.description = seq.description.replace("_R_gi", "gi")
			logger.debug("Rewriting _R_gi in desc. After change its %s" % seq.description)

	backup_filename = fasta_filename + ".bak"
	logger.debug("Backing up and restandarizing %s . Backup is %s" % (fasta_filename, backup_filename))
	shutil.copy(fasta_filename, backup_filename)
	write_standarize_fasta(fasta_filename, seqs)

def runGblocks(input_fasta,UserFlag_dict):
	#Keep OriginalAligned File:
	shutil.copy(input_fasta,input_fasta+'_original')

	num_of_seqs=0
	header_dict={} # for reconstruction
	ShortDesc_records=[]
	for seq_record in SeqIO.parse(input_fasta, "fasta"):
		header_dict[seq_record.id]=seq_record.description
		seq_record.description = seq_record.id
		ShortDesc_records.append(seq_record)
		num_of_seqs+=1
	#The file on which we'll run Gblocks (short header lines):
	t_temp=open(input_fasta+'_temp.fasta','w')
	SeqIO.write(ShortDesc_records,t_temp,"fasta")
	t_temp.close()

	#Now run Gblocks:
	gb_b1 = int(num_of_seqs/2) + 1
	gb_b2 = int(0.85 * num_of_seqs)
	gb_b3 = 8
	gb_b4 = 10
	gb_b5 = 'n'
	if UserFlag_dict['gb_SmallerBlocks'] == 'true': gb_b4 = 5
	if UserFlag_dict['gb_AllowGaps'] == 'true': gb_b5 = 'h'
	if UserFlag_dict['gb_LessStrictFlanking'] == 'true': gb_b2 = gb_b1

	Gblocks_cmd = 'Gblocks %s -t=d b1=%d -b2=%d -b3=%d -b4=%d -b5=%s' %(input_fasta+'_temp.fasta',gb_b1,gb_b2,gb_b3,gb_b4,gb_b5)
	logger.debug("Performing Gblocks command: %s" %Gblocks_cmd)
	cluster_dir = os.path.dirname(input_fasta)
	os.chdir(cluster_dir)
	#exec_external_command_redirect_output(Gblocks_cmd, cluster_dir + "/gBlocks.out", cluster_dir + "/gBlocks.err")
	os.system(Gblocks_cmd)
	output_gblocks = input_fasta+'_temp.fasta'+'-gb'  #name given by Gblocks

	#Create file with original headers, same name as input file:
	Final_records=[]
	for seq_record in SeqIO.parse(output_gblocks, "fasta"):
		seq_record.description = header_dict[seq_record.id]
		Final_records.append(seq_record)
	f_out=open(input_fasta,'w')
	SeqIO.write(Final_records,f_out,"fasta")

	return

def runTrimal(input_fasta,Trimal_cf):
	#Keep OriginalAligned File:
	shutil.copy(input_fasta,input_fasta+'_original')
	filtered_fasta = input_fasta + '_filtered'

	#Keep original headers since Triml cut them short:
	header_dict={}
	for seq_record in SeqIO.parse(input_fasta, "fasta"):
		header_dict[seq_record.id]=seq_record.description

	#Now run Trimal:
	Triml_path = ott_config['diff_soft']['TrimlPath']
	trimal_cmd = Triml_path + 'trimal -in %s -out %s -gt %s' %(input_fasta,filtered_fasta,Trimal_cf)
	logger.debug("Running cmd: %s" %trimal_cmd)
	os.system(trimal_cmd)

	#Create file with original headers, same name as input file:
	Final_records=[]
	for seq_record in SeqIO.parse(filtered_fasta, "fasta"):
		seq_record.description = header_dict[seq_record.id]
		Final_records.append(seq_record)
	f_out=open(input_fasta,'w')
	SeqIO.write(Final_records,f_out,"fasta")

	return

def align_seq_file(fasta_filename, output_aligned_fasta_filename,context):
	#parameters for Alignment:
	mafft_iter = context.UserFlags_dict['MAFFT_maxiterate']
	mafft_pairMethod = context.UserFlags_dict['PairwiseAlignmentMethod']

	if mafft_pairMethod == 'mafft_globalpair':
		mafft_cline = MafftCommandline(input=fasta_filename, nuc=True, quiet=True, maxiterate=int(mafft_iter), adjustdirection=True, globalpair=True)
	elif mafft_pairMethod == 'mafft_localpair':
		mafft_cline = MafftCommandline(input=fasta_filename, nuc=True, quiet=True, maxiterate=int(mafft_iter), adjustdirection=True, localpair=True)
	elif mafft_pairMethod == 'mafft_genafpair':
		mafft_cline = MafftCommandline(input=fasta_filename, nuc=True, quiet=True, maxiterate=int(mafft_iter), adjustdirection=True, genafpair=True)
	else:
		#Default Command
		mafft_cline = MafftCommandline(input=fasta_filename, nuc=True, adjustdirection=True, quiet=True)
	logger.info("Building alignment. Executing the following command - %s", mafft_cline)
	stdout, stderr = mafft_cline()
	handle = open(output_aligned_fasta_filename, "w")
	handle.write(stdout)
	handle.close()

	#Run Filter method to avoid large msa files:
	if context.UserFlags_dict['FilterMSA_Method'] == 'Trimal':
		logger.debug("MSA Filter method: Trimal")
		Trimal_cf = context.UserFlags_dict['Trimal_CutOff']
		runTrimal(output_aligned_fasta_filename,Trimal_cf)
		logger.info("Files after Trimal at: %s" %output_aligned_fasta_filename)
	elif context.UserFlags_dict['FilterMSA_Method'] == 'Gblocks':
		logger.debug("MSA Filter method: Gblocks")
		runGblocks(output_aligned_fasta_filename,context.UserFlags_dict)

	#sys.exit()

	#Original CODE -
	#mafft_cline = MafftCommandline(input=fasta_filename, maxiterate=1000, localpair=True)
	#logger.info("Building alignment. Executing the following command - %s", mafft_cline)
	#stdout, stderr = mafft_cline()
	#handle = open(output_aligned_fasta_filename, "w")
	#handle.write(stdout)
	#handle.close()

	#remove After debug
	#MSA_cmd = 'mafft --nuc --quiet --adjustdirection %s > %s' %(fasta_filename,output_aligned_fasta_filename)
	#exec_external_command_redirect_output(command_to_exec=MSA_cmd, outfile= local_dir + "/mafft_log.out",errfile=local_dir + "/mafft_log.err")
	#mafft_cline = MafftCommandline(input=fasta_filename, maxiterate=1000, localpair=True)
	#logger.info("Building alignment. Executing the following command - %s", mafft_cline)
	#logger.info("Building alignment. Executing the following command - %s", MSA_cmd)
	#stdout, stderr = mafft_cline()
	#handle = open(output_aligned_fasta_filename, "w")
	#handle.write(stdout)
	#handle.close()


def get_species_names_in_counts(counts_filename):
	logger.debug("Getting species names from counts  %s" % counts_filename)
	species_name_list = set()
	all_records = list(SeqIO.parse(counts_filename, "fasta"))
	for rec in all_records:
		species_name_list.add(rec.description)
	return species_name_list


def get_species_names_in_tree(trees_filename):
	logger.debug("Getting species names from trees %s" % trees_filename)
	species_name_list = set()
	trees = list(Phylo.parse(trees_filename, "newick"))
	for tree in trees:
		for clade in tree.get_terminals():
			species_name_list.add(clade.name)
	return species_name_list