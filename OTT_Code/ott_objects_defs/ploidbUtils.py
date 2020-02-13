import time
import pandas as pd
import os
import zipfile
import copy
from ete3 import Tree
from ete3 import NCBITaxa
from collections import Counter
import random
from Bio import SeqIO
from Bio import Phylo

from ott_objects_defs.ObjectJSONEncoder import ObjectJSONEncoder
from ott_objects_defs.PloiDbContext import PloiDbContext
from ott_objects_defs.ploidbCommon import *
from ott_objects_defs.handleMultAccessions import get_seq_with_max_average_blast_score

from ott_objects_defs.PathHelper import PathHelper
from ott_scripts.NR_NameMatch import do_resolve_names



__author__ = 'Moshe'

#--------------------------------------------------------------------------------------------------------
def statusFail_LogFile(context, status_line):
	logger.info(
		"=======================================================================================================")
	logger.info("===                 %s	         " % status_line)
	logger.info(
		"=======================================================================================================")
	# Result indication file:
	with open(context.final_status, 'a') as f_status:
		f_status.write("Failed - %s \n" % status_line)
	return


#--------------------------------------------------------------------------------------------------------
def get_organism_from_fasta(file_toMergeTo):
	species_list = []
	with open(file_toMergeTo, "rU") as main_f:
		for line in main_f:
			if '>' in line:
				line_rep = line.replace('>', '')
				name = line_rep.rstrip()
				species_list.append(name)
	return species_list

#--------------------------------------------------------------------------------------------------------
def is_empty_file(f_handle):

	for line in f_handle:
		if not line.strip():
			continue
		else:
			#Not empty
			return 1
	#Empty File
	return 0

#--------------------------------------------------------------------------------------------------------
def get_parent_taxId(name_binomial_val,parent_species_names_by_parent_tax_id_dict):

	for taxId, name in parent_species_names_by_parent_tax_id_dict.items():
		if name == name_binomial_val:
			return taxId
	return 'Failed'

#-------------------------------------------------------------------------------------------------------
def time_estimation_for_user(context):

	Tree_method = context.UserFlags_dict['Tree_Method']
	Ouutgroup_flag = context.UserFlags_dict['Outgroup_Flag']
	if Tree_method == 'ExaML' or Tree_method == 'RAxML':
		if Ouutgroup_flag == 'None':
			first_var = -6.07
			sec_var = -0.21
			trd_var = 0.11
		else:
			first_var = 99.34
			sec_var = 0.34
			trd_var =  0.14
	else:
		#MrBayes
		if Ouutgroup_flag == 'None':
			first_var = -14.09
			sec_var = -0.01
			trd_var = 0.2
		else:
			first_var = 273.34
			sec_var = 1.49
			trd_var = 0.27

		#y ~ 273.34 + 1.49 * num_taxa + 0.27 * num_seqs  (R2 = 0.57)

	num_taxa = len(context.species_list_names)
	num_seqs = count_seq_Infasta(context.fasta_seq_filename)

	calc = first_var + (sec_var*num_taxa)+(trd_var*num_seqs)


	with open(context.working_dir+'/calc_time_vars.txt','w') as f_time_calc:
	# if Guidance /triml/gblocks and/or Bootstrap options were selected:
		if context.UserFlags_dict['FilterMSA_Method'] != 'None' or context.UserFlags_dict['Use_BootS_RxOn'] != 'Off' or context.UserFlags_dict['Use_BootS_ExOn'] != 'Off':
			f_time_calc.write(', Currently runtime estimations are not provided when using MSA filtering tools or the bootstrap option')
		else:
			logger.debug("Calculating running time estimation:")
			logger.debug(str(first_var) + ' + ' + str(sec_var) + ' * ' +str(num_taxa) + ' + ' + str(trd_var) + ' * ' + str(num_seqs))
			if calc <= 5:
				f_time_calc.write(', estimated running time: less than 5min')
			elif calc > 100:
					hr = int(calc / 60);
					min = calc - (hr * 60)
					if min > 30:
						hr+=1
						f_time_calc.write(', estimated running time: ~%dh' %(hr))
					else:
						f_time_calc.write(', estimated running time: ~%dh' %(hr))
			else:
				f_time_calc.write(', estimated running time: ~%.0fmin' %calc)

	return
#-------------------------------------------------------------------------------------------------------
def export_dict_to_csv(dict, csv_name, headings):
	with open(csv_name, "w", newline='') as fp:
		w = csv.writer(fp, dialect='excel', quoting=csv.QUOTE_NONNUMERIC)
		w.writerow(headings)
		for key in dict.keys():
			w.writerow([key] + dict[key])

#-------------------------------------------------------------------------------------------------------
def import_csv_to_dict(csv_name):
	dict ={}
	with open(csv_name, "r") as fp:
		r = csv.reader(fp)
		for x in r:
			dict[x[0]] = x[1:]
	dict.pop("cluster_ID")
	return dict

#----------------------------------------------------------------------------------------------------
# prev used in blastSequencesLocaly as:
# from giListFromFasta import giListFromFastaFile, gi_list_str_from_fasta_file
#
# given a fasta file - write the GI to an output file
#
def giListFromFastaFile(fastaFilename, giListOutFilename):
	logger.info("about to get GI list from " + fastaFilename)
	giListOutHandle = open(giListOutFilename, 'w')
	fastaHandle = open(fastaFilename, "r")
	for seq_record in SeqIO.parse(fastaHandle, "fasta"):
		gi = getPropertyFromFastaSeqHeader(seq_record.description, "gi")
		giListOutHandle.write(gi + "\n")

	giListOutHandle.close()

# given a fasta file - return all the GIs in the fasta separated by comma
# michal: added ' around the gi since it has a structure of 'ZX123123123.2' (after replacing the gi number with the accession id)
def gi_list_str_from_fasta_file(fastaFilename):
	logger.info("about to get GI list from " + fastaFilename)
	fastaHandle = open(fastaFilename, "r")
	gis_str = ""
	for seq_record in SeqIO.parse(fastaHandle, "fasta"):
		gi = getPropertyFromFastaSeqHeader(seq_record.description, "gi")
		gis_str += ", " + "'" + gi + "'"

	return gis_str[2:]
#----------------------------------------------------------------------------------------------------

def removed_after_clustering(context):

	#read the groups file in the clustering directoryand check which taxa is not included in clusters
	#that have more than 5 species in it:


	return

#-------------------------------------------------------------------------------------------------------
def append_seq_to_largeTaxId(context):
#write_standarize_fasta(context.largeSynToAdd_file, LargeTaxon_edited_seqs)

	syn_of_large_seqs = list(SeqIO.parse(context.largeSynToAdd_file, "fasta"))

	for seq in syn_of_large_seqs:
		taxId = getPropertyFromFastaSeqHeader(seq.description, "taxonid") #this is already the large taxa
		large_fasta_file = context.largeSpeciesDir + '/seqs_for_LargeTaxId_' + str(taxId) + '.fasta'
		with open(large_fasta_file,'a') as f_toWrite:
			f_toWrite.write(">%s\n%s\n" % (seq.description, seq.seq))
			#SeqIO.write(seq, f_toWrite, "fasta")
			logger.info("Large taxon %s was added with sequences of a synonym/subsp" %taxId)

	return

#----------------------------------------------------------------------------------------------------
def check_methods_diff(context,original_run_dir):

	params_original_run = original_run_dir + '/params.txt'
	with open(params_original_run) as params_origin_f:
		for line in params_origin_f:
			logger.debug("Line of original params file: %s" %line)
			if 'MSA_Software' in line:
				msa_original_meth = line.split(':')[1]
				logger.debug("Original MSA software is: %s" % msa_original_meth)
				logger.debug("Current MSA soft is: %s" % context.UserFlags_dict['MSA_Software'])
				if msa_original_meth == context.UserFlags_dict['MSA_Software']:
					context.align_method_diff = 'false'
				else:
					context.align_method_diff = 'true'
				return

#----------------------------------------------------------------------------------------------------
def is_outgroup_included(genus_outgroup,Alignment_file):

	with open(Alignment_file,'r') as fasta_f:
		for line in fasta_f:
			if genus_outgroup in line:
				return 'Yes'
	return 'No'

#----------------------------------------------------------------------------------------------------
#This function will check if the gi str already in the cluster and then this large taxa will be skipped, since we can't add 2 spc with the same name to the concat file:
def Check_if_gi_InCluster(gi_str,f_no_multiple_accessions):

	logger.debug("Check for taxonid already in file to avoids large spc addition: %s " %gi_str)
	with open(f_no_multiple_accessions,'r') as f_fasta:
		for line in f_fasta:
			if '>' in line:
				taxId_str = '|taxonid|'+gi_str+'|'
				if taxId_str in line:
					return 'skip'
	return 'no_skip'


#----------------------------------------------------------------------------------------------------

#This function will take a list of names and check if they are scientific names on ncbi
#If so it will return thier taxa ids:
def check_scientific_on_ncbi(name_vs_acc_name_dict):

	logger.debug("Number of input names from Plant list:")
	logger.debug(len(name_vs_acc_name_dict.keys()))
	name_vs_acc_name_Scien_dict={}
	name_vs_ncbiTaxId_dict={}
	#check using ete3:
	names_without_taxId=[]
	names_with_taxId={}

	db_name = ott_config['general']['DB_DIR'] + ott_config['concat_headir']['NCBI_NAMES_NR']
	conn = sqlite3.connect(db_name)
	curs = conn.cursor()
	db_query_str = turn_list_to_DB_list(name_vs_acc_name_dict.keys())
	query_for_scientificCheck = "select TaxID,Name,Type from ncbi_names_db where lower(Name) in (%s)" % (db_query_str)
	logger.debug("query_for_scientificCheck")
	logger.debug(query_for_scientificCheck)
	curs.execute(query_for_scientificCheck)
	rows = curs.fetchall()
	for row in rows:
		TaxId = row[0]
		Name = row[1]
		Type = row[2]
		if Type == 'scientific name':
			names_with_taxId[Name] = TaxId
			name_vs_acc_name_Scien_dict[Name] = name_vs_acc_name_dict[Name]
			name_vs_ncbiTaxId_dict[Name] = TaxId
		else:
			logger.debug("NAME not scientific %s" % Name)
			names_without_taxId.append(Name)

	logger.debug("Number of names found in NCBI (as scientific)")
	logger.debug(len(names_with_taxId.keys()))
	return name_vs_acc_name_Scien_dict,name_vs_ncbiTaxId_dict

#----------------------------------------------------------------------------------------------------

def get_all_species_and_Accepted(context,db_filename,high_ranked_list,Names_AfterNR_list,id_vs_accepted_ids_dict,writer_nameTax):

#select ID from TPL_all where genus in('Crassulaceae','Abacopteris')
#union
#select ID from TPL_all where family in ('Crassulaceae','Abacopteris')
	names_list=[]
	acc_names_acc_taxid_dict = {}

	conn = sqlite3.connect(db_filename)
	curs = conn.cursor()


	#Get data for high ranked names - all species ids under the high ranked names given:
	ids_list=[]
	db_query_str = turn_list_to_DB_list(high_ranked_list)
	if context.UserFlags_dict['NR_DB_name'] == 'PLT':
		query_for_highRanked_Ids = "select ID,\"Accepted ID\" from TPL_all where lower(genus) in(%s) and \"Taxonomic status in TPL\" is not 'Unresolved' union select ID,\"Accepted ID\" from TPL_all where lower(family) in(%s) and \"Taxonomic status in TPL\" is not 'Unresolved'" % (db_query_str,db_query_str)
	elif context.UserFlags_dict['NR_DB_name'] == 'COL':
		# Check which are relevant for COL - High ranks--> kingdom phylum class order superfamily family genericName genus subgenus
		# kingdom , phylum , class , "order" , superfamily , family , genus , subgenus
		query_for_highRanked_Ids = "select taxonID,acceptedNameUsageID from col_names_table where lower(kingdom) in(%s) or lower(phylum) in(%s) or lower(class) in(%s) or lower(\"order\") in(%s) or lower(family) in(%s) or lower(genus) in(%s)" % (db_query_str,db_query_str,db_query_str,db_query_str,db_query_str,db_query_str)
	logger.debug("query_for_highRanked_Ids")
	logger.debug(query_for_highRanked_Ids)
	if db_query_str:
		curs.execute(query_for_highRanked_Ids)
		rows = curs.fetchall()
		for row in rows:
			logger.debug(row)
			id = str(row[0]).strip()
			if str(row[1]).strip() == 'None' or str(row[1]).strip() == 'nan':
				id_acc = id
			else:
				id_acc = str(row[1]).strip()
			if id not in id_vs_accepted_ids_dict.keys():
				id_vs_accepted_ids_dict[id]=id_acc
	logger.debug("Id vs Id acc after adding the high ranked species to id_vs_accepted_ids_dict")
	logger.debug(id_vs_accepted_ids_dict)

	#Get data for high ranked species ids - ids and names vs their accepted ones:
	db_query_str = turn_list_to_DB_list(id_vs_accepted_ids_dict.keys())
	if context.UserFlags_dict['NR_DB_name'] == 'PLT':
		query_hr_syn_vs_acc= "SELECT Tpl_ID,Specie,Accepted_ID,Accepted_Specie from TPL_Name_vs_Accepted where TPL_ID in (%s)" % (db_query_str)
	elif context.UserFlags_dict['NR_DB_name'] == 'COL':
		query_hr_syn_vs_acc= "SELECT tmp1.taxonID, tmp1.scientificName, tmp1.acceptedNameUsageID, tmp2.AccName FROM \
				(  SELECT taxonID,scientificName,acceptedNameUsageID from col_names_table where taxonID in (%s)  )tmp1 \
				LEFT JOIN \
				( SELECT taxonID,scientificName as AccName from col_names_table where taxonID in (SELECT acceptedNameUsageID from col_names_table where taxonID in (%s)) )tmp2 \
				ON tmp1.acceptedNameUsageID = tmp2.taxonID" % (db_query_str,db_query_str)
	id_vs_acc_id_dict={}
	name_vs_acc_name_dict={}
	logger.debug("query_hr_syn_vs_acc")
	logger.debug(query_hr_syn_vs_acc)
	if db_query_str:
		curs.execute(query_hr_syn_vs_acc)
		rows = curs.fetchall()
		for line in rows:
			logger.debug(line)
			Tpl_ID=line[0]
			Specie=line[1].rstrip()
			if line[2]:
				Accepted_ID=line[2]
			else:
				Accepted_ID=Tpl_ID
			if line[3]:
				Accepted_Specie=line[3].rstrip()
			else:
				Accepted_Specie=Specie
			names_list.append(Tpl_ID)
			names_list.append(Specie)
			id_vs_acc_id_dict[Tpl_ID]=Accepted_ID
			name_vs_acc_name_dict[Specie]=Accepted_Specie
	logger.debug("All Date ids, names, ids dict and names dict:")
	logger.debug(names_list)
	logger.debug(id_vs_acc_id_dict)
	logger.debug(name_vs_acc_name_dict)


	##Get data for species names - ids and names vs their accepted ones:
	#db_query_str = turn_list_to_DB_list(Names_AfterNR_list)
	#query_syn_vs_acc= "SELECT Tpl_ID,Specie,Accepted_ID,Accepted_Specie from TPL_Name_vs_Accepted where lower(Specie) in (%s)" % (db_query_str)
	#logger.debug("query_syn_vs_acc")
	#logger.debug(query_syn_vs_acc)
	#curs.execute(query_syn_vs_acc)
	#rows = curs.fetchall()
	#for line in rows:
	#	Tpl_ID=line[0]
	#	Specie=line[1]
	#	Accepted_ID=line[2]
	#	Accepted_Specie=line[3]
	#	id_vs_acc_id_dict[Tpl_ID]=Accepted_ID
	#	name_vs_acc_name_dict[Specie]=Accepted_Specie
	#	names_list.append(Specie)

	#Now check how many names are scientific names on NCBI:
	###name_vs_ncbiScien_TaxId_dict = check_scientific_on_ncbi(name_vs_acc_name_dict)
	name_vs_acc_Scientific_dict, name_vs_ncbiScien_TaxId_dict = check_scientific_on_ncbi(name_vs_acc_name_dict)

	logger.debug("name_vs_acc_Scientific_dict")
	logger.debug(name_vs_acc_Scientific_dict)

	#Check all names in the list: if the accepted name has a scientific name on NCBI and therfor has a TaxID, get this taxID,
	#Otherwise, generate a fiction taxID for this accepted name for the name replacement
	for name in name_vs_acc_Scientific_dict.keys():
		acc_name = name_vs_acc_Scientific_dict[name]
		if acc_name not in name_vs_ncbiScien_TaxId_dict.keys(): #couldn't find TaxId
			#Check if Accepted name is already in the dictionary so we wouldn't create a 2nd random id:
			if acc_name in acc_names_acc_taxid_dict.keys():
				fiction_taxId = acc_names_acc_taxid_dict[acc_name]
			else:
				fiction_taxId = int('999'+str(random.randint(0,10))+str(random.randint(0,10))+str(random.randint(0,10))+str(random.randint(0,10))+str(random.randint(0,10))+str(random.randint(0,10))+'999')
				acc_names_acc_taxid_dict[acc_name]=str(fiction_taxId)
				#Check if this will solev the error when Matched PLT - (example : Lithospermum)
				context.species_names_by_tax_id[str(fiction_taxId)] = acc_name
			logger.debug("Create fictional Accepted TaxId for plant list without a scientific NCBI name: %s, %s" %(acc_name,fiction_taxId))
		else:
			acc_names_acc_taxid_dict[acc_name] = name_vs_ncbiScien_TaxId_dict[acc_name]

	logger.debug("Final dictionary: Name Vs Accpeted name, length %d" %len(name_vs_acc_Scientific_dict.keys()))
	logger.debug(name_vs_acc_Scientific_dict)	#names Vs Accpeted names
	logger.debug("Final dictionary: Accepted names Vs TaxIds, length %d" %len(acc_names_acc_taxid_dict.keys()))
	logger.debug(acc_names_acc_taxid_dict)		#Accepted names TaxIds
	logger.debug("Accepted names with real TaxIds, length %d" %len(name_vs_ncbiScien_TaxId_dict.keys()))
	logger.debug(name_vs_ncbiScien_TaxId_dict)

	for key in name_vs_ncbiScien_TaxId_dict.keys():
		logger.debug("%s,%s" %(key,name_vs_ncbiScien_TaxId_dict[key]))

	#Update dictionaries for Merging Syn with Accepted names:
	for name in name_vs_acc_Scientific_dict.keys():
		Acc_name = name_vs_acc_Scientific_dict[name]
		TaxId = name_vs_ncbiScien_TaxId_dict[name]
		Acc_TaxID = acc_names_acc_taxid_dict[Acc_name]
		#update OTT dict and list:
		context.syn_acc_TaxID_dict[str(TaxId)] = str(Acc_TaxID)
		context.syn_acc_dict[str(Acc_TaxID)] = Acc_name
		context.species_list_ids.append(str(TaxId))
		context.species_list_names.append(name)
		context.species_names_by_tax_id[str(TaxId)]=name
		writer_nameTax.writerow([str(TaxId), name])


	logger.debug("context.syn_acc_TaxID_dict   Syn Vs Accepted PLT")
	logger.debug(context.syn_acc_TaxID_dict)		#Accepted names TaxIds
	logger.debug("syn_acc_dict PLT")
	logger.debug(context.syn_acc_dict)		#Accepted names TaxIds

	conn.close()
	logger.debug("END of DEBUG section !!!!!!!!!!!")
	#sys.exit()
	return

#Name res/Matched names for Plant list (in future implement on other DBs as well:
def resolve_names_plt(context,db_name,Names_input_list):

	names_taxid_f=open(context.species_names_taxid_file,mode="w",encoding='utf8',newline='')
	writer_name_tax = csv.writer(names_taxid_f,delimiter=',')

	#In case the list includes sigle names - >High ranked, we need to check if family or genus and retrive all relevant taxa:
	high_ranked_list=[]
	names_for_nr=[]
	species_ids=[]
	id_vs_acc_id_dict=[]
	name_vs_acc_name_dict=[]
	for name in Names_input_list:
		if ' ' not in name:
			high_ranked_list.append(name)
		else:
			names_for_nr.append(name)
	Names_AfterNR_list,id_vs_accepted_ids_dict = name_plt_matching(db_name,names_for_nr,context)
	logger.debug('Names_AfterNR_list')
	logger.debug(Names_AfterNR_list)
	logger.debug('id_vs_accepted_ids_dict')
	logger.debug(id_vs_accepted_ids_dict)
	get_all_species_and_Accepted(context,db_name,high_ranked_list,Names_AfterNR_list,id_vs_accepted_ids_dict,writer_name_tax)

	#Set the lists for OTT pipeline:
	logger.debug("NEW NEW NEW:")
	logger.debug(context.species_list_ids)
	logger.debug(context.species_list_names)
	logger.debug(context.species_names_by_tax_id)

	#Once we have the TaxID list we need to extract all subsp and their parents for merging:
	logger.debug("SYNONYM:")
	logger.debug(context.syn_acc_TaxID_dict)
	logger.debug(context.syn_acc_dict)


	return


# Usign ete to get all descendants - NCBI (including better handling of the intra-spesific group:
def list_added_syn_acc(context,db_filename,accepted_ids_list):
	all_syn_list=[]
	Syn_Acc_dict={}
	db_query_str = turn_list_to_DB_list(accepted_ids_list)
	conn = sqlite3.connect(db_filename)
	curs = conn.cursor()
	if 'COL' in db_filename:
		query_for_all_syn = "SELECT scientificName from col_names_table where acceptedNameUsageID IN (%s) or taxonID IN (%s)" % (db_query_str,db_query_str)
	else:
		query_for_all_syn = "SELECT Specie from TPL_Name_vs_Accepted where Accepted_ID IN (%s) or Tpl_ID IN (%s)" % (db_query_str,db_query_str)
	logger.debug("db_query_str")
	logger.debug(db_query_str)
	logger.debug(query_for_all_syn)
	curs.execute(query_for_all_syn)
	rows = curs.fetchall()
	with open(context.working_dir+'/Syn_names.txt','w') as f_syn:
		for row in rows:
			name=list(row)[0].rstrip()
			logger.debug(name)
			f_syn.write(name+'\n')
			all_syn_list.append(name)
	logger.debug(len(rows))
	#Get Accepted names for all names:
	names_list_query=turn_list_to_DB_list(all_syn_list)
	if 'COL' in db_filename:
		query_for_all_acc= "SELECT tmp1.scientificName, tmp2.AccName FROM \
				(  SELECT taxonID,scientificName,acceptedNameUsageID from col_names_table where taxonID in (%s)  )tmp1 \
				LEFT JOIN \
				( SELECT taxonID,scientificName as AccName from col_names_table where taxonID in (SELECT acceptedNameUsageID from col_names_table where taxonID in (%s)) )tmp2 \
				ON tmp1.acceptedNameUsageID = tmp2.taxonID" % (names_list_query,names_list_query)
	else:
		query_for_all_acc= "SELECT Specie,Accepted_Specie from TPL_Name_vs_Accepted where lower(Specie) IN (%s)" % (names_list_query)
	logger.debug(query_for_all_acc)
	curs.execute(query_for_all_acc)
	rows = curs.fetchall()
	with open(context.working_dir+'/Acc_names.txt','w') as f_acc:
		f_acc.write("Name,Acc_Name\n")
		for row in rows:
			name=list(row)[0].rstrip()
			ac_name=list(row)[1].rstrip()
			f_acc.write(name+','+ac_name+'\n')
			logger.debug(name+','+ac_name+'\n')
			Syn_Acc_dict[name]=ac_name

	conn.close()

	return all_syn_list,Syn_Acc_dict

#----------------------------------------------------------------------------------------
def break_ncbi_line(line): #return Spc_Rank,Spc_TaxID,Spc_ParentTaxID,Spc_Name,Spc_Type
	Spc_Rank = str(line[0])
	Spc_TaxID = str(line[1])
	Spc_ParentTaxID = str(line[2])
	Spc_Name = line[3].replace(",","").replace("'","").replace("-"," ").replace("[","").replace("]","")
	Spc_Type = str(line[5])
	return Spc_Rank,Spc_TaxID,Spc_ParentTaxID,Spc_Name,Spc_Type

def update_parents_names_dict(context,db_name):

	parents_taxIds_list=[]
	for taxId in context.parent_species_tax_id_by_subsp_tax_id.values():
		parents_taxIds_list.append(taxId)

	TaxID_parents_list_str = turn_list_to_DB_list(parents_taxIds_list)
	query = "select * from ncbi_names_db where TaxID IN (%s) and Type is 'scientific name'" %TaxID_parents_list_str
	db_rows = query_ncbi_db(query,db_name)
	for line in db_rows:
		Spc_Rank,Spc_TaxID,Spc_ParentTaxID,Spc_Name,Spc_Type = break_ncbi_line(line)
		context.parent_species_names_by_parent_tax_id[Spc_TaxID] = Spc_Name
	return



def read_csv_to_dict(csv_file):
	# open the file in universal line ending mode
	with open(csv_file, 'rU') as infile:
		# read the file as a dictionary for each row ({header : value})
		reader = csv.DictReader(infile)
		data = {}
		for row in reader:
			for header, value in row.items():
				try:
					data[header].append(value)
				except KeyError:
					data[header] = [value]
	return data



def compare_taxa_files(original_final_species_f,new_species_input_f):

	update_species_list=[];remove_species_list=[]
	list_origin = [];list_new = []
	with open (original_final_species_f,'r') as f_old:
		for line in f_old:
			line= line.strip()
			list_origin.append(line)
			logger.debug(line)
	with open (new_species_input_f,'r') as f_new:
		for line in f_new:
			line= line.strip()
			list_new.append(line)
			logger.debug(line)

	#Get the common species - to keep in the new list:
	common_species_list = list(set(list_origin).intersection(list_new))
	#List of species to add:
	add_species_list = [a for a in list_new if (a not in list_origin)]    # Check which species to add
	remove_species_list = [a for a in list_origin if (a not in list_new)]    # Check which species to remove

	#logger.debug("---------------------   Rerun Edit taxa List data    -------------------------\n")

	update_species_list = common_species_list + add_species_list
	logger.debug("-------> Updated New List: ")
	logger.debug(update_species_list)
	logger.debug("------------------------------------------------------------------------------\n")
	#return:  update_species_list,remove_species_list,addNew_species_list
	return update_species_list,remove_species_list,add_species_list


def zipdir(path, zipHandle):
	for root, dirs, files in os.walk(path):
		for file in files:
			ziph.write(os.path.join(root, file))


def merge_seq_files(seq_files_to_merge, merged_seqs_filename,its_flag):

	if its_flag == 'ITS':
		logger.info("merging ITS the following files %s into %s " % (seq_files_to_merge[1], merged_seqs_filename))
		addFrag_cmd = 'mafft --addfragments ' + seq_files_to_merge[1] + ' --multipair ' + merged_seqs_filename + ' > ' + merged_seqs_filename + '_added'
		os.system(addFrag_cmd)
		shutil.copyfile(merged_seqs_filename + '_added',merged_seqs_filename)
		return
	else:
		logger.info("merging the following files %s into %s " % (','.join(seq_files_to_merge), merged_seqs_filename))
		seq_list = list()
		for f in seq_files_to_merge:
			for record in SeqIO.parse(f, "fasta"):
				seq_list.append(record)
		merge_seq_handle = open(merged_seqs_filename, 'w')
		SeqIO.write(seq_list, merge_seq_handle, "fasta")
		merge_seq_handle.close()
		return

def remove_seq_from_file(specie_toRemove, fasta_file_toEdit): #227904673
	logger.info("Remove '%s' from %s " % (specie_toRemove, fasta_file_toEdit))
	seq_list = list()
	for record in SeqIO.parse(fasta_file_toEdit, "fasta"):
		if specie_toRemove in record.description:
			logger.debug("Removing %s, seq: %s " %(specie_toRemove,record.description))
		else:
			seq_list.append(record)
	edit_seq_handle = open(fasta_file_toEdit, 'w')
	SeqIO.write(seq_list, edit_seq_handle, "fasta")
	edit_seq_handle.close()


def get_pickled_filename(subdir):
	l = list()
	for filename in os.listdir(subdir):
		#logger.debug("next file %s" % filename)
		if filename.endswith(PloiDbContext.PICKLED_FILE_EXTENSION):
			l.append(os.path.join(subdir,filename))
	if len(l) > 1:
		print("warn - more than one pickled file found in %s returning %s" % (subdir,l[0]))
		return l[0]
	elif len(l) == 0:
		print("Error - no pickled file found in %s " % subdir)
		return None
	else:
		return l[0]


def pickle_context(context,pickled_filename,pickled_json_filename):
	with open(pickled_filename, "wb") as handle:
		pickle.dump(context,handle)
		logger.info("Created pickled file %s out of the context with %d clusters" % (pickled_filename,len(context.cluster_contexts)))
	with open(pickled_json_filename, "w") as handle:
		json.dump(context,handle, sort_keys = True, indent = 4,cls=ObjectJSONEncoder)
		logger.info("Created pickled file %s out of the context with %d clusters" % (pickled_filename,len(context.cluster_contexts)))


#---------------------------------------------------------------------------------------------------------------------------
# Polling Jmodel jobs to check if they are done:
def call_Check_Jmt_Done(jobs_dirs_list):

	sleep_period = 2*60	# 30min
	# Check status of jobs in: jmt_JobsList.txt every 10min:
	jmt_polling_time = ott_config['general']['JMT_PollingTime']
	jmt_start_time = time.time()
	end_time = jmt_start_time + (int(jmt_polling_time) * 3600)
	while time.time() < end_time:
		dirs_count=len(jobs_dirs_list); jmt_done_count=0
		logger.info("Polling JMT...time -> %s" % time.time())
		#Check if done:
		for dir_path in jobs_dirs_list:
			logger.info("Checking dir %s" % dir_path)
			if os.path.exists(dir_path+'/JMT_DONE'):
				logger.info("Found JMT_DONE at %s" % dir_path)
				jmt_done_count+=1
		if dirs_count == jmt_done_count:
			return 0
		else:
			time.sleep(sleep_period)
	return 1


#---------------------------------------------------------------------------------------------------------------------------
# Polling Bootstrap runs in parallel:
def polling_BS_Done(XML_Dir,bootstrap_val):

	Flag_done='True'

	sleep_period = 10*60	# check every 10 min
	# Check status of jobs in: jmt_JobsList.txt every 10min:
	jmt_polling_time = 120
	jmt_start_time = time.time()
	end_time = jmt_start_time + (int(jmt_polling_time) * 3600)
	while time.time() < end_time:
		done_cnt=0
		for Boots_idx in range(0, int(bootstrap_val)):
			BS_idx_Dir=XML_Dir+'/BS_'+str(Boots_idx)
			if os.path.exists(BS_idx_Dir+'/ExaML_modelFile.examlOut'):
				logger.debug("File: %s Exists" %(BS_idx_Dir+'/ExaML_modelFile.examlOut'))
				done_cnt+=1
		logger.info("Polling BS...time -> %s, done for this point %s out of %s" % (time.time(),done_cnt,str(bootstrap_val)))
		#Check if done:
		if done_cnt==int(bootstrap_val):
			return 0
		else:
			time.sleep(sleep_period)
	return 1
#---------------------------------------------------------------------------------------------------------------------------

def check_cluster_BlastAll(cluster_num,seq_id,working_dir):
	blast_header = ['queryId', 'subjectId', 'percIdentity', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 'subjectEnd', 'eVal', 'bitScore']
	path_blast_allVall = working_dir + "/" + str(cluster_num) + "/blast_all-v-all_total.blastn"
	# Read blast all file into a data-structure:
	blast_all_DataStruct = pd.read_csv(path_blast_allVall, sep='\t',names = blast_header)
	# filter only the relevant results for selected seq_id:
	#filtered = blast_all_DataStruct[(blast_all_DataStruct['queryId'] == seq_id)]
	if 'seqid' in seq_id:
		seq_id_filt = seq_id.split('seqid|')[1]
	else:
		seq_id_filt=seq_id
	filtered = blast_all_DataStruct[(blast_all_DataStruct['queryId'] == seq_id_filt)]
	if not filtered.empty:
		# Drop duplicates of subjectId -> take only the first data (the high score):
		No_duple = filtered.drop_duplicates(subset=['subjectId'], keep='first')
		if not No_duple.empty:
			# 10% cutoff value:
			location_of_cutoff = int(0.1*len(No_duple))
			CuttOff_val = sorted(No_duple['bitScore'].values)[location_of_cutoff]
			logger.debug("check_cluster_BlastAll: seqid-%s CutOff-%s" %(seq_id,CuttOff_val))
			return CuttOff_val
		else:
			logger.debug("check_cluster_BlastAll: Empty No_duple after Drop duplicates of subjectId, seq_id %s" %seq_id)
			return 999
	else:
		logger.debug("check_cluster_BlastAll: Empty structure after filter of relevant results for selected seq_id %s" %seq_id)
		return 999

#---------------------------------------------------------------------------------------------------------------------------
#Used for creating pruned tree that will only contain species with CCDB counts data
def get_species_ToRemove(genus,working_dir):

	Species_withCounts_list=[]; Species_withCounts_count = 0; species_name='Null'
	species_list_OnTree_list=[];species_list_OnTree_count=0
	Species_ToRemove_list=[];Species_ToRemove_count=0
	#Open file to write species names for removal:
	species_remove_list_file = working_dir + 'concat/remove_species_list.txt'
	f_species_to_remove = open(species_remove_list_file,'w')

	# Get all species names with counts (according to counts file based on CCDB):
	counts_file = COUNTS_DIR + genus + '.counts'
	f= open(counts_file,'r')
	for line in f:
		line=line.strip('\n')
		if '>' in line:
			species_name = line.strip('>')
			Species_withCounts_list.append(species_name)

	# Check which species on Tree, by searching the concat file:
	concat_file = working_dir + 'concat/' + genus + '-concat-aligned.fasta'
	f_concat= open(concat_file,'r')
	for line in f_concat:
		line=line.strip('\n')
		if '>' in line:
			species_name = line.strip('>')
			species_list_OnTree_list.append(species_name)
			if species_name not in Species_withCounts_list:
				#List of species to remove: on tree but no counts
				Species_ToRemove_list.append(species_name)
				f_species_to_remove.write(species_name+'\n')
				Species_ToRemove_count+=1

	species_list_OnTree_count = len(species_list_OnTree_list)

	# Check if qualifies to run chromevol: more than 5 species with counts:
	if (species_list_OnTree_count - Species_ToRemove_count) < 5:
		f_lessThanFive = open(working_dir + 'concat/LessThanFiveWithCounts.txt','w')
		f_lessThanFive.close()

	return

#---------------------------------------------------------------------------------------------------------------------------
def replaceToOrganismName(f_origin,organism_name):
	# with is like your try .. finally block in this case
	with open(f_origin, 'r') as file:
		# read a list of lines into data
		data = file.readlines()

	# now change the 2nd line, note that you have to add a newline
	data[0] = '>' + organism_name +'\n'

	# and write everything back
	with open(f_origin, 'w') as file:
		file.writelines( data )
#---------------------------------------------------------------------------------------------------------------------------
def replaceToLowerCase(f_origin):
	# with is like your try .. finally block in this case
	with open(f_origin, 'r') as file:
		# read a list of lines into data
		data = file.readlines()

	# now change the 2nd line, note that you have to add a newline
	max_idx=len(data)
	idx_num = 0
	while idx_num < max_idx:
		if '>' in data[idx_num]:
			step='step'
			idx_num+=1
		else:
			data[idx_num] = data[idx_num].lower()
			idx_num+=1

	# and write everything back
	with open(f_origin, 'w') as file:
		file.writelines( data )

#---------------------------------------------------------------------------------------------------------------------------

def remove_shortLong_seqs(blast_results, db_fasta_file, path_of_selected_seq, blast_FilterOutput,log_file):
	#remove_shortLong_seqs(blast_results, db_fasta_file, path_of_selected_seq, blast_FilterOutput,log_file)
	script_path = ott_config['general']['OTT_MAIN'] + 'ott_scripts/'
	# Preprocessing the blast results in order to filter sequences that are not covered by a certain threshold in the blast results
	logger.debug("Remove Short/Long sequences from blast results:")
	os.system ("python " + script_path + "FilterBlast_LengthParam.py -i " + blast_results + " -i_seq " + path_of_selected_seq + " -f " + db_fasta_file + " -o " + blast_FilterOutput +" > " + log_file)

#---------------------------------------------------------------------------------------------------------------------------
def checkLargeSeqBlast(out_blastn):
	# send back the first rresult after filtering results of Too long/short sequences:
	with open(out_blastn,'r') as f_blastAll:
		first_line = f_blastAll.readline()
		line_split = first_line.split('\t')
		#list_results = [line_split[1],format(float(line_split[2]),'.2f'),format(float(line_split[3]),'.2f'),format(float(line_split[10]),'.2f'),format(float(line_split[11]),'.2f')]
		list_results = [line_split[1],format(float(line_split[11]),'.2f')]
		return list_results

#---------------------------------------------------------------------------------------------------------------------------
def add_seq_to_AlignedFile(new_seq,Aligned_file):

	# Example command :   "mafft --add new_sequences --reorder existing_alignment > output"
	Aligned_file_out = Aligned_file + '_added'
	mafft_add_cmd = "mafft --add " + new_seq + " --reorder " + Aligned_file + " > " + Aligned_file_out
	os.system(mafft_add_cmd)
	#replace names so the new fasta will have the original name:
	Aligned_file_origin = Aligned_file + '_origin'
	os.rename(Aligned_file, Aligned_file_origin)
	os.rename(Aligned_file_out, Aligned_file)

	#os.rename(Aligned_file_out, Aligned_file)

# where this is true trees and models are not created - this is mainly to verify the entire pipeline is OK without having to wait for tree creation
# simulation_mode = False
#
#---------------------------------------------------------------------------------------------------------------------------
def Check_and_Add_LargeTaxIds(context):

	largeSpeciesTaxIDsList=[]
	f_fastaLargeList = open(context.largeSpeciesFileList,'r')
	for line in f_fastaLargeList:
		start_TaxId_str='LargeTaxId_'
		end_TaxId_str='.fasta'
		TaxId_num = line[line.find(start_TaxId_str)+len(start_TaxId_str):line.rfind(end_TaxId_str)]
		largeSpeciesTaxIDsList.append(TaxId_num)

	if len(largeSpeciesTaxIDsList) > 0:
		logger.debug("TaxIds with many sequences on Genbank - get sequences and compare to existing clusters:")
		f_added_TaxId_gis = open(context.summary_dir + '/Added_LargeTaxIds.txt', 'w')
	for LargeTaxId in largeSpeciesTaxIDsList:
		logger.debug("Handeling Large TaxId %s" %str(LargeTaxId))
		file_of_LargeSeqs = context.largeSpeciesDir + '/seqs_for_LargeTaxId_' + str(LargeTaxId) + '.fasta'
		#1. Check Large Species and create db using formatdb:
		# Create data base for comparison of selected sequences in each cluster:
		# if db exists don't run formatdb:
		if not os.path.exists(file_of_LargeSeqs+'.nhr'):
			formatdb_cmd = "formatdb -i " + file_of_LargeSeqs + " -pF"  #formatdb -i %s -pF
			os.system(formatdb_cmd)
		Dict_Add_gi_to_Cluster={} #keys are clusters, values are gi's to add

		with open(context.working_dir + "/LargeSpecies/LOG_blastall_" + str(LargeTaxId) + ".txt",'w') as log_file:
			log_file.write("ClusterNum,BitScore_cutOff,BestBitScore_NewSpecie,NewSpecie_gi\n")
			for cluster_context in context.cluster_contexts:
				illegal_seqs=[]
				#2.
				i_clust = cluster_context.index
				if os.path.exists(context.working_dir + "/" + str(i_clust) + '/ITS_CLUSTER'):
					logger.debug("ITS cluster, skip Large species (done separately)")
					continue
				if os.path.exists(context.working_dir + "/" + str(i_clust) + '/ITS_CLUSTER_2'):
					logger.debug("ITS2 cluster, skip Large species (done separately)")
					continue
				if os.path.exists(context.working_dir + "/" + str(i_clust) + '/ITS_CLUSTER_3'):
					logger.debug("ITS3 cluster, skip Large species (done separately)")
					continue
				#for i_clust in num_of_clusters:
				path_of_selected_seq = context.working_dir + "/" + str(i_clust) + "/seq-to-blast"
				with open(path_of_selected_seq,'r') as f_seq:
					logger.debug("Get data of selected seq of cluster %s" %str(i_clust))
					first_line_seq = f_seq.readline()
					start_seqid='seqid|'
					end_seqid='|description|'
					seq_id = first_line_seq[first_line_seq.find(start_seqid)+len(start_seqid):first_line_seq.rfind(end_seqid)]
					start_gi = '>gi|'
					end_gi = '|taxonid|'
					gi_num = first_line_seq[first_line_seq.find(start_gi)+len(start_gi):first_line_seq.rfind(end_gi)]
					logger.debug("Selected seq: gi %s, seqid %s" %(gi_num,seq_id))

				blast_results = context.working_dir + "/LargeSpecies/" +str(LargeTaxId) + "_seq_vs_LargeSeq_" + str(i_clust) + ".blastn"
				blast_cmd = "blastall -p blastn -d " + file_of_LargeSeqs + " -i " + path_of_selected_seq + " -v 150000 -b 150000 -e 0.1 -m 8 > " + blast_results   # -m 9 will give a header to the output file
				logger.debug("Running blast cmd: %s" % blast_cmd)
				os.system(blast_cmd)
				#Filter results:  blast_results, db_fasta_file, blast_FilterOutput,log_file
				db_fasta_file = file_of_LargeSeqs
				#Filter blast results of Too Long/Short sequences:
				#cluster_context.blast_results_filename
				blast_FilterOutput = context.working_dir + "/LargeSpecies/seq_vs_LargeSeqFiltered_" + str(i_clust) + ".blastn"
				log_file_filter = context.working_dir + "/LargeSpecies/log_filter_"+ str(i_clust)+".txt"
				remove_shortLong_seqs(blast_results, db_fasta_file, path_of_selected_seq, blast_FilterOutput,log_file_filter)
				cutOff_BitScore = check_cluster_BlastAll(i_clust,seq_id,context.working_dir)
				if os.stat(blast_FilterOutput).st_size == 0 or cutOff_BitScore == 999:
					largeSeq_blastResults = [0,0]
				else:
					largeSeq_blastResults = checkLargeSeqBlast(blast_FilterOutput)
				#Add seq to aligned file:
				log_file.write(str(i_clust) + ',' + str(cutOff_BitScore)+ ',' + str(largeSeq_blastResults[1]) + ',' + str(largeSeq_blastResults[0])+'\n')
				#Check if BitScore is >= to CutOff value:
				if (float(largeSeq_blastResults[1]) >= cutOff_BitScore):
					f_added_TaxId_gis.write("\nTaxID-%s: " % str(LargeTaxId))
					if largeSeq_blastResults[0] in Dict_Add_gi_to_Cluster.values():
						illegal_seqs.append(largeSeq_blastResults[0])
					else:
						Dict_Add_gi_to_Cluster[i_clust] = largeSeq_blastResults[0]
						f_added_TaxId_gis.write("%s - %s," %(str(i_clust),str(largeSeq_blastResults[0])))
						gi_num = largeSeq_blastResults[0].split('|')[1]
						#Create seq file for adding:
						add_seq_file = context.working_dir + "/LargeSpecies/add_taxId_" + str(LargeTaxId) + "_gi_" + str(gi_num) + "_toClust_"+ str(i_clust) + '.fasta'
						add_seq_withDesc_file = context.working_dir + "/LargeSpecies/add_desc_taxId_" + str(LargeTaxId) + "_gi_" + str(gi_num) + "_toClust_"+ str(i_clust) + '.fasta'
						with open(file_of_LargeSeqs, "rU") as fasta_large_f:
							all_seqs = list(SeqIO.parse(fasta_large_f, "fasta"))
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
					fasta_large_f.close()
			#Check each cluster was selected once:
			logger.debug("Dict_Add_gi_to_Cluster:")
			logger.debug(Dict_Add_gi_to_Cluster)
			if illegal_seqs:
				logger.debug("Large TaxID matched more than one cluster, need to verify different gi's !!")
			else:
				for cluster_context in context.cluster_contexts:
					i_clust = cluster_context.index
					if i_clust in Dict_Add_gi_to_Cluster.keys():
						gi_num = Dict_Add_gi_to_Cluster[i_clust].split('|')[1]
						#Check if cluster already includes this gi, if so skip it:
						Skip_flag = Check_if_gi_InCluster(str(LargeTaxId),cluster_context.all_seqs_fasta_filename_no_multiple_accessions)
						logger.debug("Skip_flag is %s" %Skip_flag)
						if Skip_flag != 'skip':
							add_seq_file = context.largeSpeciesDir + "/add_taxId_" + str(LargeTaxId) + "_gi_" + str(gi_num) + "_toClust_"+ str(i_clust) + '.fasta'
							add_seq_withDesc_file = context.working_dir + "/LargeSpecies/add_desc_taxId_" + str(LargeTaxId) + "_gi_" + str(gi_num) + "_toClust_"+ str(i_clust) + '.fasta'
							#Add seq to existing all seqs file:
							shutil.copyfile(cluster_context.all_seqs_fasta_filename_no_multiple_accessions,cluster_context.all_seqs_fasta_filename_no_multiple_accessions+'_original')
							its_file_name = context.working_dir + '/' + str(i_clust) + '/ITS_CLUSTER'
							if os.path.exists(its_file_name):
								logger.debug("Cluster %s is an ITS cluster")
							#if cluster_context.is_its_cluster:
								#merge to already aligned file:
								merge_seq_files([cluster_context.all_seqs_fasta_filename_no_multiple_accessions,
											 add_seq_withDesc_file],
											cluster_context.all_seqs_fasta_filename_no_multiple_accessions,'ITS')
							else:
								merge_seq_files([cluster_context.all_seqs_fasta_filename_no_multiple_accessions,
											 add_seq_withDesc_file],
											cluster_context.all_seqs_fasta_filename_no_multiple_accessions,'Not_ITS')
							replaceToLowerCase(cluster_context.all_seqs_fasta_filename_no_multiple_accessions)
						else:
							logger.debug("skipped large species on cluster %s since taxon %s already in it" %(str(i_clust),gi_num))
		logger.debug("Illegal seqs:")
		logger.debug(illegal_seqs)

#------------------------------------------------------------------------------------------------------------------
def remove_its_seqs(context,fasta_newTax_filename,fasta_newTax_ITS_f,fasta_newTax_WithoutITS_f):

	# Filtering the non ITS sequences
	scripts_path = context.cluster_script_dir + "/clustering/"
	output_dir = context.rerun_new_taxIds_dir

	#Split new taxID into 2 files: without ITS and ITS_only:
	remove_none_its_command = "perl " + scripts_path + "getSeqsWithITS.pl " + fasta_newTax_filename + " " + fasta_newTax_ITS_f
	logger.info("Calling getSeqsWithITS (for rerun new species ITS clustering) - " + remove_none_its_command)
	retval = exec_external_command_redirect_output(command_to_exec=remove_none_its_command,
												   outfile=output_dir + "/ITS_only/getSeqsWithITS.out",
												   errfile=output_dir + "/ITS_only/getSeqsWithITS.err")

	remove_its_command = "perl " + scripts_path + "getSeqsWithoutITS.pl " + fasta_newTax_filename + " " + fasta_newTax_WithoutITS_f
	logger.info("Calling getSeqsWithoutITS (for ITS removal) - " + remove_its_command)

	exec_external_command_redirect_output(command_to_exec=remove_its_command,
										  outfile=output_dir + "/Without_ITS/getSeqsWithoutITS.out",
										  errfile=output_dir + "/Without_ITS/getSeqsWithoutITS.err")

	# Check if files are empty and if so delete them:
	if os.stat(fasta_newTax_ITS_f).st_size == 0:
		os.remove(fasta_newTax_ITS_f)
	else:
		context.rerun_new_its_files_list.append(fasta_newTax_ITS_f)
		formatdb_cmd = "formatdb -i " + fasta_newTax_ITS_f + " -pF -o T"
		os.system(formatdb_cmd)
	if os.stat(fasta_newTax_WithoutITS_f).st_size == 0:
		os.remove(fasta_newTax_WithoutITS_f)
	else:
		context.rerun_new_non_its_files_list.append(fasta_newTax_WithoutITS_f)
		formatdb_cmd = "formatdb -i " + fasta_newTax_WithoutITS_f + " -pF -o T"
		os.system(formatdb_cmd)

	return
#------------------------------------------------------------------------------------------------------------------
def Add_ITS_seqs(context):
	logger.debug("Add_ITS_seqs() is EMPTY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ")
	# In case the original run has a combined ITS file, we can run the usual way:
	combined_original = context.UserFlags_dict['OriginJobDir'] + '/cluster_its/combined.fasta'
	if os.stat(combined_original).st_size == 0:
		return False
	else:
		return True


#-----------------------------------------------------------------------------------------------------------
def rerun_get_representative_seq(taxon_id,taxon_fasta_filename,taxon_selected_filename,taxon_blast_filename):

	logger.info("taxonid = %s - handling multiple accessions" % taxon_id)
	representative_seq = get_seq_with_max_average_blast_score(taxon_fasta_filename,taxon_blast_filename)

	with open(taxon_selected_filename, "w") as handle:
		SeqIO.write(representative_seq, handle, "fasta")





#------------------------------------------------------------------------------------------------------------------------
#2017 New Name resolve function for OneTwoTree:
# This function should inable the merge option of sub-species with theie parents and synonyms with their accepted oraganism
#------------------------------------------------------------------------------------------------------------------------
def resolve_names(context):

	# Get all taxIDs from fasta file:
	logger.info("Resolving names")
	all_seqs = list(SeqIO.parse(context.fasta_seq_filename, "fasta"))

	ids_to_resolved_names = dict()
	ids_not_resolved = list()

	resolved_seqs = list()
	resolved_names_taxonid = dict()  # since we need that the taxonid
	unresolved_names = set()

	for seq in all_seqs:
		# logger.debug("Processing seq.id=%s seq.description=%s" % (seq.id, seq.description))
		original_taxonid = getPropertyFromFastaSeqHeader(seq.description, "taxonid")
		original_name = getPropertyFromFastaSeqHeader(seq.description, "organism")

		skip_seq = False
		resolved_organism_name = context.species_names_by_tax_id[original_taxonid]
		if resolved_organism_name not in resolved_names_taxonid:
			resolved_names_taxonid[resolved_organism_name] = original_taxonid

		#MD - Update species list to include synonyms - For outgroup selection reasons:
		#if original_name not in context.species_list_names:
		#	context.species_list_names.append(original_name)
		restandarize_seq(seq, original_organism=original_name, organism=resolved_organism_name,taxonid=resolved_names_taxonid[resolved_organism_name],
						 original_taxonid=original_taxonid)  # setting seq properties
		resolved_seqs.append(seq)


	shutil.copyfile(context.fasta_seq_filename, context.fasta_seq_org_names_filename)
	write_standarize_fasta(context.fasta_seq_filename, resolved_seqs)
	context.unresolved_species_names = unresolved_names
	logger.debug("All unresolved species are: ")
	logger.debug(unresolved_names)
	logger.debug("All species including synonyms: ")
	logger.debug(context.species_list_names)



def Create_NR_input(context,taxa_list):

	InputFile_NR = open(context.working_dir+'/InputFor_NR.csv', 'w')
	InputFile_NR.write("Id,Name\n")
	index=1
	for taxa in taxa_list:
		InputFile_NR.write('%s,%s\n' %(index,taxa))
		index+=1
	return


#--------------------------------------------------def-------------------------------------------------------
def query_ncbi_db(query,db_name):

	db_name = ott_config['general']['DB_DIR'] + ott_config['concat_headir']['NCBI_NAMES_NR']
	logger.debug("Perform the following DB query: %s" %query)

	#connect to DB
	conn_ncbi = sqlite3.connect(db_name)
	curs_ncbi = conn_ncbi.cursor()
	curs_ncbi.execute(query)
	rows_ncbi = curs_ncbi.fetchall()
	return rows_ncbi

#--------------------------------------------------def-------------------------------------------------------
def extract_names_from_rows(context,rows_db,db_name):

	TaxID_List=[]
	synTaxID_List=[]
	Names_List=[]
	synNames_List=[]
	TaxId_rank_dict=dict()
	Name_type_dict=dict()
	TaxId_name_dict=dict()
	TaxId_Parent_dict=dict()
	#['Rank','TaxID','ParentTaxID','Name','Other','Type'])
	for line in rows_db:
		Rank=line[0]
		TaxID=line[1]
		ParentTaxID=line[2]
		Name=line[3].replace(",","").replace("'","").replace("-"," ").replace("[","").replace("]","")
		Other=line[4]
		Type=line[5]
		#In case of Name match/Resolution, split scientific and synonyms:
		if context.UserFlags_dict['NameResType'] == 'SpellingCorrection' or context.UserFlags_dict['NameResType'] == 'MatchedNames':
			if Type == 'scientific name':
				logger.debug(' ')
			#Keep synonym data - all other name types will be classified as syn names for this specific scientific name:
			else:#
				synTaxID_List.append(str(TaxID))
				synNames_List.append(Name)
				context.synTaxId_rank_dict[str(TaxID)]=Rank
				#context.synName_type_dict[Name]=Type
				context.synTaxId_name_dict[str(TaxID)]=Name
				context.synTaxId_Parent_dict[str(TaxID)]=str(ParentTaxID)
		TaxID_List.append(str(TaxID))
		Names_List.append(Name)
		context.TaxId_rank_dict[str(TaxID)]=Rank
		context.Name_type_dict[Name]=Type
		context.TaxId_name_dict[str(TaxID)]=Name
		context.Name_TaxId_dict[Name]=str(TaxID)
		context.TaxId_Parent_dict[str(TaxID)]=str(ParentTaxID)

	#print(TaxID_List)
	TaxId_rank_dict = context.TaxId_rank_dict.copy()
	Name_type_dict = context.Name_type_dict.copy()
	TaxId_name_dict = context.TaxId_name_dict.copy()
	Name_TaxId_dict = context.Name_TaxId_dict.copy()
	TaxId_Parent_dict = context.TaxId_Parent_dict.copy()

	return TaxID_List,Names_List,TaxId_rank_dict,Name_type_dict,TaxId_name_dict,Name_TaxId_dict,TaxId_Parent_dict

#--------------------------------------------------def-------------------------------------------------------
def extract_names_from_rows_debug(context,rows_db,db_name):

	TaxID_List=[]
	synTaxID_List=[]
	Names_List=[]
	synNames_List=[]
	TaxId_rank_dict=dict()
	TaxId_name_dict=dict()
	Name_TaxId_dict=dict()
	TaxId_Parent_dict=dict()
	#['Rank','TaxID','ParentTaxID','Name','Other','Type'])
	for line in rows_db:
		Rank=line[0]
		TaxID=line[1]
		ParentTaxID=line[2]
		Name=line[3].replace(",","").replace("'","").replace("[","").replace("]","")
		Other=line[4]
		Type=line[5]
		#In case of Name match/Resolution, split scientific and synonyms:
		if context.UserFlags_dict['NameResType'] == 'SpellingCorrection' or context.UserFlags_dict['NameResType'] == 'MatchedNames':
			if Type == 'scientific name':
				logger.debug(' ')
			#Keep synonym data - all other name types will be classified as syn names for this specific scientific name:
			else:#
				synTaxID_List.append(str(TaxID))
				synNames_List.append(Name)
				context.synTaxId_rank_dict[str(TaxID)]=Rank
				context.synTaxId_name_dict[str(TaxID)]=Name
				context.synTaxId_Parent_dict[str(TaxID)]=str(ParentTaxID)
		TaxID_List.append(str(TaxID))
		Names_List.append(Name)
		context.TaxId_rank_dict[str(TaxID)]=Rank
		context.TaxId_name_dict[str(TaxID)]=Name
		context.Name_TaxId_dict[Name]=str(TaxID)
		context.TaxId_Parent_dict[str(TaxID)]=str(ParentTaxID)

	#print(TaxID_List)
	TaxId_rank_dict = context.TaxId_rank_dict.copy()
	TaxId_name_dict = context.TaxId_name_dict.copy()
	Name_TaxId_dict = context.Name_TaxId_dict.copy()
	TaxId_Parent_dict = context.TaxId_Parent_dict.copy()

	return TaxID_List,Names_List,TaxId_rank_dict,TaxId_name_dict,Name_TaxId_dict,TaxId_Parent_dict

#--------------------------------------------------def-------------------------------------------------------
def turn_list_to_DB_list(list_name):

	db_names_str=''
	for item in list_name:
		name_lowercase = str(item).lower()
		item_add="'"+name_lowercase+"',"
		db_names_str+=item_add
	final_str = db_names_str[:-1]
	return final_str
#--------------------------------------------------def-------------------------------------------------------
def Split_low_high_Ranks(context,input_TaxId_list,db_name):

	High_rank_TaxId=[]
	Spc_rank_TaxId=[]
	TaxId_list_query=turn_list_to_DB_list(input_TaxId_list)
	#print(TaxId_list_query)
	#query = "select * from ncbi_names_db where ParentTaxID IN (%s) and Type is 'scientific name'" %TaxId_list_query
	query = "select * from ncbi_names_db where ParentTaxID IN (%s) and Type not in ('common name','authority')" %TaxId_list_query
	db_rows = query_ncbi_db(query,db_name)
	TaxId_list,Names_list,TaxId_rank_dict,Name_type_dict,TaxID_Name_dict,TaxId_Parent_dict = extract_names_from_rows(context,db_rows,db_name)
	for TaxId in TaxId_list:
		if context.TaxId_rank_dict[TaxId] in ['species','subspecies','varietas']:
			Spc_rank_TaxId.append(TaxId)
		else:
			if TaxId not in High_rank_TaxId:
				High_rank_TaxId.append(TaxId)   #To check in case species have kids with rank 'sub' or such
	#Check for TaxIds that are not Parents -> Already species level:
	#TaxId_withoutKids = list(set(input_TaxId_list) - set(TaxId_list))
	#logger.debug("Species who have no Kids:")
	#logger.debug(TaxId_withoutKids)
	#for TaxId_spc in TaxId_withoutKids:
	#	Spc_rank_TaxId.append(TaxId_spc)

	return Spc_rank_TaxId,High_rank_TaxId
#--------------------------------------------------def-------------------------------------------------------
def Split_low_high_Ranks_debug(context,input_TaxId_list,db_name):

	High_rank_TaxId=[]
	Spc_rank_TaxId=[]
	TaxId_list_query=turn_list_to_DB_list(input_TaxId_list)
	#print(TaxId_list_query)
	#query = "select * from ncbi_names_db where ParentTaxID IN (%s) and Type is 'scientific name'" %TaxId_list_query
	query = "select * from ncbi_names_db where ParentTaxID IN (%s)" %TaxId_list_query
	db_rows = query_ncbi_db(query,db_name)
	TaxId_list,Names_list,TaxId_rank_dict,TaxID_Name_dict,Name_TaxId_dict,TaxId_Parent_dict = extract_names_from_rows_debug(context,db_rows,db_name)
	for TaxId in TaxId_list:
		if context.TaxId_rank_dict[TaxId] in ['species','subspecies','varietas']:
			Spc_rank_TaxId.append(TaxId)
		#else:  # cancelled this else so species who have kids will also be included in the parents list
		if TaxId not in High_rank_TaxId:
			High_rank_TaxId.append(TaxId)   #To check in case species have kids with rank 'sub' or such
	#Check for TaxIds that are not Parents -> Already species level:
	TaxId_withoutKids = list(set(input_TaxId_list) - set(TaxId_list))
	logger.debug("Species who have no Kids:")
	logger.debug(TaxId_withoutKids)
	for TaxId_spc in TaxId_withoutKids:
		Spc_rank_TaxId.append(TaxId_spc)

	return Spc_rank_TaxId,High_rank_TaxId
#--------------------------------------------------def-------------------------------------------------------
def create_file_for_NR(Names_input_list,context):

	f_nr = open(context.working_dir+'/NR_input.csv','w')
	f_nr.write("Id,Name\n")
	idx=1
	for Name in Names_input_list:
		if Name:
			f_nr.write("%d,%s\n" %(idx,Name))
			idx+=1
	f_nr.close()
	return


#--------------------------------------------------def-------------------------------------------------------
# Check for matched names - Spelling correction only:
def db_name_matching(db_name,Names_input_list,context):
	create_file_for_NR(Names_input_list,context)
	nr_match_cmd = "python " + context.ott_scripts_dir + "/NR_NameMatch.py \
	--db-filename %s \
	--input-filename %s \
	--results-filename %s \
	--log-filename %s \
	--authfield False" %(db_name,context.working_dir+'/NR_input.csv',context.working_dir+'/NR_output.csv',context.working_dir+'/NR_Log.txt')

	nr_log_out = context.working_dir+'/NR_LOG.OE'
	nr_log_err = context.working_dir+'/NR_LOG.ER'
	exec_external_command_redirect_output(nr_match_cmd, nr_log_out, nr_log_err)


	matched_names_dict = {}
	matched_names_scores_dict = {}
	original_names_list = []
	no_match_names = []
	new_list=[]
	accepted_ids_list = []
	# Read Name res output file and extract the matched names:
	reader = csv.reader(open(context.working_dir+'/NR_output.csv', 'r'))
	# Add file to see for each name what happened ith it: no match at all, poor match(<0.8) etc...
	NR_analyze = open(context.nr_analyze,'w')
	NR_analyze.write('OriginalName,MatchedName,Score,includeFlag\n')#if IncludeFlag is 0 than this name is not included in our pipeline
	for i, rows in enumerate(reader):
		if i == 0: continue
	#for rows in reader:
		id = rows[0]
		parent_id = rows[1]
		Type = rows[2]
		if Type == 'Synonym' or Type == 'Unresolved' or Type == 'synonym':
			if parent_id != '-':
				accepted_ids_list.append(parent_id)
		elif Type == 'Accepted' or Type == 'accepted name':
			if id != '-':
				accepted_ids_list.append(id)
		matched_name = rows[10]
		score = float(rows[9])
		original_name = rows[5]
		original_names_list.append(original_name)
		matched_names_dict[original_name] = matched_name
		matched_names_scores_dict[original_name] = score
		if float(matched_names_scores_dict[original_name]) >= 0.8:
			new_list.append(matched_name)
			NR_analyze.write('%s,%s,%.2f,1\n' %(original_name,matched_name,score))
		else:
			#Names that will be excluded from this run
			no_match_names.append(original_name)
			logger.debug(original_name)
			logger.debug(matched_names_scores_dict[original_name])
			if float(matched_names_scores_dict[original_name]) == 0:
				NR_analyze.write('%s,%s,%.2f,0\n' %(original_name,matched_name,score))

	logger.debug("NR statistics:")
	logger.debug("Input list included: %d species" %(i))
	logger.debug("Final list include: %d species" %len(new_list))
	logger.debug("Removed species (low match score): %d species" %len(no_match_names))

	# Original Input vs NR output: search for names that were not found in the NR (no match)
	not_in_nr = list(set(Names_input_list)-set(original_names_list))
	no_match_names.extend(not_in_nr)

	logger.debug("MISSING NAMES:")
	logger.debug(no_match_names)
	logger.debug(matched_names_dict)
	logger.debug(matched_names_scores_dict)

	#New: in case the user is using another DB (not NCBI) we need to handle the synonym names differently. Currently only plant list:
	#All syn and Accepted names of a syn/acc name in the input list according to the selected DB will be extracted and added to the input list:
	if context.UserFlags_dict['NR_DB_name'] == 'PLT':
		db_file = ott_config['general']['DB_DIR'] + ott_config['concat_headir']['TPL_ALL_DB']
		list_with_syn_list,syn_vs_acc_dict = list_added_syn_acc(context,db_file,accepted_ids_list)
		logger.debug(list_with_syn_list)
		logger.debug(syn_vs_acc_dict)
		return list_with_syn_list,syn_vs_acc_dict
	elif context.UserFlags_dict['NR_DB_name'] == 'COL':
		db_file = ott_config['general']['DB_DIR'] + ott_config['concat_headir']['COL_NAMES_DB']
		list_with_syn_list,syn_vs_acc_dict = list_added_syn_acc(context,db_file,accepted_ids_list)
		logger.debug(list_with_syn_list)
		logger.debug(syn_vs_acc_dict)
		return list_with_syn_list,syn_vs_acc_dict
	#Create Acc/Syn lists:

	logger.debug("new_list:")
	logger.debug(new_list)

	return new_list
	#return list_with_syn_list,syn_vs_acc_dict


#--------------------------------------------------def-------------------------------------------------------
# Check for matched names - Spelling correction only:
def name_plt_matching(db_name,Names_input_list,context):
	create_file_for_NR(Names_input_list,context)
	nr_match_cmd = "python " + context.ott_scripts_dir + "/NR_NameMatch.py \
	--db-filename %s \
	--input-filename %s \
	--results-filename %s \
	--log-filename %s \
	--authfield False" %(db_name,context.working_dir+'/NR_input.csv',context.working_dir+'/NR_output.csv',context.working_dir+'/NR_Log.txt')

	nr_log_out = context.working_dir+'/NR_LOG.OE'
	nr_log_err = context.working_dir+'/NR_LOG.ER'
	exec_external_command_redirect_output(nr_match_cmd, nr_log_out, nr_log_err)


	matched_names_dict = {}
	matched_names_scores_dict = {}
	original_names_list = []
	no_match_names = []
	unresolved_names=[]
	new_list=[]
	accepted_ids_list = []
	id_vs_accepted_ids_dict = {}
	# Read Name res output file and extract the matched names:
	reader = csv.reader(open(context.working_dir+'/NR_output.csv', 'r'))
	# Add file to see for each name what happened ith it: no match at all, poor match(<0.8) etc...
	NR_analyze = open(context.nr_analyze,'w')
	NR_analyze.write('OriginalName,MatchedName,Score,includeFlag\n')#if IncludeFlag is 0 than this name is not included in our pipeline
	for i, rows in enumerate(reader):
		if i == 0: continue
	#for rows in reader:
		id = rows[0]
		parent_id = rows[1]
		Type = rows[2]
		matched_name = rows[10]
		score = float(rows[9])
		original_name = rows[5]
		#names_after_nr = []
		original_names_list.append(original_name)
		matched_names_dict[original_name] = matched_name
		matched_names_scores_dict[original_name] = score
		# only if the name is valid (according to score) save accepted data and matched name:
		if float(matched_names_scores_dict[original_name]) >= 0.8:
			if Type == 'Synonym' or Type == 'synonym':
				new_list.append(matched_name)
				NR_analyze.write('%s,%s,%.2f,1,Synonym\n' %(original_name,matched_name,score))
				if parent_id != '-':
					#accepted_ids_list.append(parent_id)
					id_vs_accepted_ids_dict[id]=parent_id
			elif Type == 'Accepted' or Type == 'accepted name' or '-': #TPL or COL
				new_list.append(matched_name)
				NR_analyze.write('%s,%s,%.2f,1,Accepted\n' %(original_name,matched_name,score))
				if id != '-':
					#accepted_ids_list.append(id)
					id_vs_accepted_ids_dict[id]=id
			elif Type == 'Unresolved' or 'Misapplied':	#unresolved names will be saved both in an unresolved list and the unmatched list
				unresolved_names.append(original_name)
				no_match_names.append(original_name)
				NR_analyze.write('%s,%s,%.2f,0,Unresolved\n' %(original_name,matched_name,score))
		else:
			no_match_names.append(original_name)
			NR_analyze.write('%s,%s,%.2f,0\n' %(original_name,matched_name,score))
			logger.debug(original_name)
			logger.debug(matched_names_scores_dict[original_name])

	logger.debug("NR statistics:")
	logger.debug("Input list included: %d species" %(i))
	logger.debug("Final list include: %d species" %len(new_list))
	logger.debug("Removed species (low match score): %d species" %len(no_match_names))
	logger.debug("Unresolved Names: %d species" %len(unresolved_names))


	# Original Input vs NR output: search for names that were not found in the NR (no match)
	not_in_nr = list(set(Names_input_list)-set(original_names_list))
	no_match_names.extend(not_in_nr)

	logger.debug("MISSING NAMES:")
	logger.debug(no_match_names)
	logger.debug(matched_names_dict)
	logger.debug(matched_names_scores_dict)

	#New: in case the user is using another DB (not NCBI) we need to handle the synonym names differently. Currently only plant list:
	#All syn and Accepted names of a syn/acc name in the input list according to the selected DB will be extracted and added to the input list:
	if context.UserFlags_dict['NR_DB_name'] == 'PLT':
		db_file = ott_config['general']['DB_DIR'] + ott_config['concat_headir']['TPL_ALL_DB']
		list_with_syn_list,syn_vs_acc_dict = list_added_syn_acc(context,db_file,accepted_ids_list)
		logger.debug(list_with_syn_list)
		logger.debug(syn_vs_acc_dict)

	#Create Acc/Syn lists:

	#return new_list
	return new_list,id_vs_accepted_ids_dict



def get_splitSpecies(split_input):

	input_split=split_input.split(',')
	Name1 = input_split[0].strip()
	Name2 = input_split[1].strip()

	return Name1,Name2


def return_msa_length(report_file):

	max_point = 0
	with open(report_file,'r') as f_report:
		firstline = True
		for line in f_report:
			if firstline == True:
				firstline=False
				continue
			columns = line.split('\t')
			if int(columns[2]) > max_point:
				max_point = int(columns[2])
	return max_point


def create_NameDate_config(context,msa_length,spl_under_name1,spl_under_name2, MinAge,MaxAge):

	origin_tree = context.summary_dir + '/Result_Tree_' + context.id + '_NotUltra.tre'
	out_tree = context.summary_dir + '/Result_Tree_' + context.id + '.tre'

	with open(context.tree_xml_namedate_config,'w') as f_config:
		f_config.write("treefile = %s\n" %origin_tree)
		f_config.write("smooth = 100\n")
		f_config.write("numsites = %d\n" %msa_length)
		f_config.write("mrca = NODE_LBL %s %s\n" %(spl_under_name1,spl_under_name2))
		f_config.write("min = NODE_LBL %s\n" %MinAge)
		f_config.write("max = NODE_LBL %s\n" %MaxAge)
		f_config.write("outfile = %s\n" %out_tree)
	f_config.close()


	return

#def create_FileNameDate_config(context,msa_length,spl_under_name1,spl_under_name2, MinAge,MaxAge):
def create_FileNameDate_config(context,msa_length,NumOfNodeDates,Name1_list,Name2_list,MinAge_list,MaxAge_list):

	origin_tree = context.summary_dir + '/Result_Tree_' + context.id + '_NotUltra.tre'
	out_tree = context.summary_dir + '/Result_Tree_' + context.id + '.tre'

	with open(context.tree_xml_namedate_config,'w') as f_config:
		f_config.write("treefile = %s\n" %origin_tree)
		f_config.write("smooth = 100\n")
		f_config.write("numsites = %d\n" %msa_length)
		node_idx=0
		while node_idx < NumOfNodeDates:
			Node_lbl =  'NODE_LBL' + str(node_idx)
			f_config.write("mrca = %s %s %s\n" %(Node_lbl,Name1_list[node_idx],Name2_list[node_idx]))
			f_config.write("min = %s %s\n" %(Node_lbl,MinAge_list[node_idx]))
			f_config.write("max = %s %s\n" %(Node_lbl,MaxAge_list[node_idx]))
			node_idx+=1
		f_config.write("outfile = %s\n" %out_tree)
	f_config.close()


	return

#--------------------------------------------------------------------------------------------------
#Checking the new way for species list reconstruction:
def return_species_list_debug(input_list,context):
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
		db_name = ott_config['general']['DB_DIR'] + ott_config['concat_headir']['NCBI_NAMES_NR']
	elif context.UserFlags_dict['NR_DB_name'] == 'PLT': # the plant list
		db_name = ott_config['general']['DB_DIR'] + ott_config['concat_headir']['TPL_ALL_DB']
	elif context.UserFlags_dict['NR_DB_name'] == 'COL': # the plant list
		db_name = ott_config['general']['DB_DIR'] + ott_config['concat_headir']['COL_NAMES_DB']
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


	if Names_input_list:
		logger.debug("User Input contains names:")
		# 2. NAME RESOLUTION section - TBD TBD !!!!!!!!!!!!!!!!!
		#----------------------------------------------------------
		#Handle the names according to user request, only relevant for species names (more than 1 str):
		#Case of Mathced names and Plant list - diff alg
		if context.UserFlags_dict['NameResType'] == 'MatchedNames' and (context.UserFlags_dict['NR_DB_name'] == 'PLT' or context.UserFlags_dict['NR_DB_name'] == 'COL'):
			#names_resolved_plt = resolve_names_plt(Names_input_list)
			resolve_names_plt(context,db_name,Names_input_list)
			return
		elif context.UserFlags_dict['NameResType'] == 'SpellingCorrection' or context.UserFlags_dict['NameResType'] == 'MatchedNames':
			for name in Names_input_list:
				if ' ' in name:
					names_for_nr.append(name)
			logger.debug("Names_input_list:")
			logger.debug(Names_input_list)
			logger.debug("names_for_nr:")
			logger.debug(names_for_nr)
			Names_input_list = list(set(Names_input_list) - set(names_for_nr))
			logger.debug("Names_input_list:")
			logger.debug(Names_input_list)
			if context.UserFlags_dict['NR_DB_name'] == 'NCBI':
				Names_AfterNR_list = db_name_matching(db_name,names_for_nr,context)
				logger.debug("Names_AfterNR_list:")
				logger.debug(Names_AfterNR_list)
				Names_input_list.extend(Names_AfterNR_list)
			else: #TPL and COL
				Names_AfterNR_list,syn_vs_acc_dict = db_name_matching(db_name,names_for_nr,context)
				Names_input_list.extend(Names_AfterNR_list)

		logger.debug("Names_input_list - After Name resolution stage:")
		logger.debug(Names_input_list)
		logger.debug("Names_AfterNR_list:")
		logger.debug(Names_AfterNR_list)

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

		logger.debug("Names_input_list")
		logger.debug(Names_input_list)
		logger.debug("context.species_names_by_tax_id")
		logger.debug(context.species_names_by_tax_id)

	# 3. Convert Names to TaxIDs:
	NamesSTRING_for_query = turn_list_to_DB_list(Names_input_list)
	query_for_Names = "select * from ncbi_names_db where lower(Name) IN (%s)" %NamesSTRING_for_query
	db_rows = query_ncbi_db(query_for_Names,db_name)
	TaxId_list,Names_list,TaxId_rank_dict,Name_type_dict,TaxID_Name_dict,Name_TaxId_dict,TaxId_Parent_dict = extract_names_from_rows(context,db_rows,db_name)
	if TaxID_input_list:
		Request_TaxIds_list = set(TaxID_input_list + TaxId_list)
	else:
		Request_TaxIds_list=list(set(TaxId_list))
	logger.debug("User requested TaxIds list:")
	logger.debug(Request_TaxIds_list)


	#combine Name HR with TaxId HR - Working with Higher ranks
	#by distinguishing between parents with spc kids and parents with higher ranked kids:
	TaxID_input_list=[]

	if Spc_Descendants == 'No':
		TaxId_list_query=turn_list_to_DB_list(Request_TaxIds_list)
		#print(TaxId_list_query)
		#query = "select * from ncbi_names_db where ParentTaxID IN (%s) and Type is 'scientific name'" %TaxId_list_query
		query = "select * from ncbi_names_db where TaxID IN (%s) and Type is 'scientific name'" %TaxId_list_query
		db_rows = query_ncbi_db(query,db_name)
		TaxID_input_list,Names_list,TaxId_rank_dict,TaxID_Name_dict,Name_TaxId_dict,TaxId_Parent_dict = extract_names_from_rows_debug(context,db_rows,db_name)
		#for item in TaxID_Name_dict.keys():
		#	logger.debug(item)
		#	logger.debug(TaxID_Name_dict[item])

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
	TaxID_input_list_str = turn_list_to_DB_list((TaxID_input_list))
	query = "select * from ncbi_names_db where TaxID IN (%s) and Type is 'scientific name'" %TaxID_input_list_str
	db_rows = query_ncbi_db(query,db_name)
	TaxId_list,Names_list,TaxId_rank_dict,Name_type_dict,TaxID_Name_dict,Name_TaxId_dict,TaxId_Parent_dict = extract_names_from_rows(context,db_rows,db_name)
	with open(context.working_dir+'/input_list_ncbi_data.csv','w') as input_data_ncbi:
		for line in db_rows:
			Spc_Rank,Spc_TaxID,Spc_ParentTaxID,Spc_Name,Spc_Type = break_ncbi_line(line)
			input_data_ncbi.write('"'+Spc_Rank+'",'+'"'+Spc_TaxID+'",'+'"'+Spc_ParentTaxID+'",'+'"'+Spc_Name+'",'+'"'+Spc_Type+'"\n')
			if Spc_Rank == 'subspecies' or Spc_Rank == 'varietas':
				context.parent_species_tax_id_by_subsp_tax_id[Spc_TaxID] = Spc_ParentTaxID
			context.Pipe_input_Spc_Names_vs_TaxId[Spc_Name] = Spc_TaxID
	#Get names of all parents for name matching:
	update_parents_names_dict(context,db_name)
	logger.debug("TaxId list - all")
	logger.debug(TaxID_Name_dict.keys())

	#Get all Accepted(for matching names) and Parents for subsp.:
	#for name in TaxID_Name_dict.values():


	#Filter final list - go over by TaxId:
	taxIds_to_Filter = list()
	logger.debug('number of species before filter - %s' %len(TaxId_rank_dict.keys()))
	if context.UserFlags_dict['Filter_SubSp'] == 'on' or context.UserFlags_dict['Filter_Hybrids'] == 'on' or context.UserFlags_dict['Filter_Unresolved'] == 'on':
		logger.debug("Perform Filter - at least one filter flag is set to 'on'")
		#Go over all TaxIds in the list:
		for taxId in TaxId_rank_dict.keys():
			name = TaxID_Name_dict[taxId]
			rank = TaxId_rank_dict[taxId]
			if context.UserFlags_dict['Filter_SubSp'] == 'on':
				if rank == 'varietas' or rank == 'subspecies' or rank == 'no rank': #Added no rank to enable all desc
					taxIds_to_Filter.append(taxId)
					logger.debug("TaxId %s will be filtered: rank is %s (intraspecific variants)" %(taxId,rank))
			if context.UserFlags_dict['Filter_Hybrids'] == 'on':
				if ' x ' in name or name.startswith("x "):	# Added Genus hybryds filter as well (begin wuth x )
					if taxId not in taxIds_to_Filter: taxIds_to_Filter.append(taxId)
					logger.debug("TaxId %s will be filtered: recognized as Hybrid, name-%s" %(taxId,name))
			if context.UserFlags_dict['Filter_Unresolved'] == 'on':
				if ' cf. ' in name or ' sp.' in name or ' aff. ' in name or ' f. ' in name:
					if taxId not in taxIds_to_Filter: taxIds_to_Filter.append(taxId)
					logger.debug("TaxId %s will be filtered: recognized as Open nomenclature, name-%s" %(taxId,name))


	#remove all Filtered TaxIds fromthe final dictionaries:
	for k in taxIds_to_Filter:
		try:
			del TaxID_Name_dict[k]
			del TaxId_rank_dict[k]
		except KeyError:
			logger.debug("Failed to remove the TaxId item %s from Names and Ranks Dictionaries" %k)
			pass

	#For debug Only:
	with open(context.working_dir+'/outNR.csv','w') as f_out:
		for key in TaxId_rank_dict.keys():
			#f_out.write('%s,%s,%s\n' %(str(key),str(TaxID_Name_dict[key]).encode("latin-1"),str(TaxId_rank_dict[key])))
			f_out.write('%s,%s,%s\n' %(str(key),str(TaxID_Name_dict[key]),str(TaxId_rank_dict[key])))
			#context.species_names_by_tax_id[str(TaxID_Name_dict[key])] = str(key)
			context.species_names_by_tax_id[str(key)] = str(TaxID_Name_dict[key])
			writer.writerow([str(key), str(TaxID_Name_dict[key])])
	f_out.close()

	#Create Dictionary to hole names for TaxId after NR:
	context.Species_names_vs_TaxId_AfterNR = copy.deepcopy(TaxID_Name_dict)

	logger.debug(context.species_names_by_tax_id)

	#Set the lists for OTT pipeline:
	context.species_list_ids = list(TaxID_Name_dict.keys())
	context.species_list_names = list(TaxID_Name_dict.values())
	#context.parent_species_tax_id_by_subsp_tax_id
	#context.parent_species_names_by_subsp_tax_id
	logger.debug("NEW NEW NEW:")
	logger.debug(context.species_list_ids)
	logger.debug(context.species_list_names)
	logger.debug(context.species_names_by_tax_id)
	#Once we have the TaxID list we need to extract all subsp and their parents for merging:
	logger.debug("SYNONYM:")
	logger.debug(context.syn_acc_TaxID_dict)
	logger.debug(context.syn_acc_dict)


	#Filter Results:
	#temp#_dict = copy.deepcopy(context.Pipe_input_Spc_Names_vs_TaxId)
	#for name in temp_dict.keys():
	#	if context.UserFlags_dict['Filter_SubSp'] == 'on':
	#		logger.debug("Filter_SubSp is on")
	#		if ' subsp ' in name or ' var ' in name:
	#			#context.Pipe_input_Spc_Names_vs_TaxId.pop(name, None)
	#			del context.Pipe_input_Spc_Names_vs_TaxId[name]
	#	if context.UserFlags_dict['Filter_Hybrids'] == 'on':
	#		logger.debug(name)
	#		if ' x ' in name:
	#			del context.Pipe_input_Spc_Names_vs_TaxId[name]
	##since we have names with both 'x' and 'cf' we need to recopy the dictionary to exclude the names already removed.
	#temp_dict = copy.deepcopy(context.Pipe_input_Spc_Names_vs_TaxId)
	#for name in temp_dict.keys():
	#	if context.UserFlags_dict['Filter_Unresolved'] == 'on':
	#		if ' cf. ' in name or ' sp. ' in name or ' aff. ' in name or ' f. ' in name:
	#			del context.Pipe_input_Spc_Names_vs_TaxId[name]
#
	#logger.debug("After Filter: Pipe_input_Spc_Names_vs_TaxId")
	#logger.debug(context.Pipe_input_Spc_Names_vs_TaxId)

	#sys.exit()

	return


#Rhis function will convert the output of BlustClust to Groups.txt file so
# we can use the orthomcl final stage of creating the clustering:

def convert_BlastClust_output_toGroups(context,PrecentIdentity,CoverageLength,blastclust_input):

	blastclust_output = context.working_dir+'/clustering/blustClust.out'
	blast_clust_cmd = 'blastclust -i %s -o %s -p F -L %s -b T -S %s -e F' %(blastclust_input,blastclust_output,CoverageLength,PrecentIdentity)
	os.system(blast_clust_cmd)
	logger.debug(blast_clust_cmd)

	f_groups_BlastClust = open(context.working_dir+'/clustering/BC_groups.txt','w')
	with open (blastclust_output,'r') as blastClust_f:
		for line in blastClust_f:
			line_split=line.strip().split(' ')
			for item in line_split: # example of an item: gi|GQ404407.1|taxonid|2903|organism|Emiliania
				acc_data = item.split('|')
				TaxId = acc_data[3]
				AcceNum = acc_data[1]
				Tax_and_Acc = str(TaxId) + '|' + str(AcceNum)
				f_groups_BlastClust.write(Tax_and_Acc + ' ')
			f_groups_BlastClust.write('\n')

	return
#-----------------------------------------------------------------------------------------------------------------
def return_leafEndsNames(tree_file_path):
	# process input from cmd
	logger.debug('returns the leaftmost and rightmost leafs')

	t = Tree(tree_file_path)
	postorder_leafs = [leaf.name for leaf in t]
	logger.debug("leaftmost leaf: %s, rightmost leaf: %s  " %(postorder_leafs[0],postorder_leafs[-1]))

	return(postorder_leafs[0],postorder_leafs[-1])
#-----------------------------------------------------------------------------------------------------------------
def create_mb_blk_constraint(context):

	os.chdir(context.concat_workdir)
	cnst_tree_f = context.working_dir+'/ConstraintTree_user.txt'
	R_Path = ott_config['diff_soft']['R_PATH']
	scripts_dir = ott_config['general']['OTT_MAIN'] + 'ott_scripts/'
	ConstraintBlk_Rcmd = (R_Path + " CMD BATCH '--args dir=" + '"' + context.concat_workdir + '"' + " tree=" + '"' + cnst_tree_f + '"' + "' " + scripts_dir + "createMrbayesConst.R")
	logger.info("Constraint R script - for MB: %s" % ConstraintBlk_Rcmd)
	os.system(ConstraintBlk_Rcmd)
	return


def	check_for_constriant_species(context):

	const_names_list=[]
	t = Tree(context.working_dir+'/ConstraintTree_user.txt', format=1)
	list_species = t.get_leaves()
	#convert to names list only:
	for item in list_species:
		logger.debug(item)
		name = str(item)
		name_cln = name.replace('\n--','').replace('_',' ')
		logger.debug(name_cln)
		const_names_list.append(name_cln)
	logger.debug("list_species from tree:")
	logger.debug(const_names_list)
	logger.debug("context.species_list_names:")
	logger.debug(context.species_list_names)

	#if list(diff.elements()):
	if set(const_names_list) == set(context.species_list_names):
		return 'Yes'
	else:
		logger.debug("Lists are not identical")
		return 'No'

#-------------------------------------------------------------------------------------
def count_taxa_in_fasta(fasta_file):

	taxa_list=[]
	with open(fasta_file,'r') as f_fasta:
		for line in f_fasta:
			if '>' in line:
				exp = r"taxonid\|(.*?)(\||\\n)"
				m = re.search(exp, line)
				if m:
					taxonId = m.group(1)
					if taxonId not in taxa_list:
						taxa_list.append(taxonId)
	num_taxa = len(set(taxa_list))

	return num_taxa,taxa_list

#-------------------------------------------------------------------------------------
def count_seq_Infasta(fasta_file):

	counter=0
	with open(fasta_file,'r') as f_fasta:
		for line in f_fasta:
			if '>' in line:
				counter+=1

	return counter

#-----------------------------------------------------------------------------------------------------------------------
def py_ITS_splitITS(TaxID,seq_desc):

	return

#
#	if 'ITS1'.lower in seq_desc.lower or 'ITS2'.lower in seq_desc.lower:
#


# given a fasta file, separate the sequences to 3 groups: (1) contain both ITS1 and ITS2, (2) contain ITS1 without ITS2, (3) contain ITS2 without ITS1
###sub splitITS{
###
###	if(@_ < 5){
###		die "usage:
###		(1) input fasta file
###		(2) output file contains all the sequences contain ITS1 and not ITS2
###		(3) output file contains all the sequences contain ITS2 and not ITS1
###		(4) output file contains all the sequences contain both ITS1 and ITS2
###		(5) genBank index
###		";
###	}
###
###	my ($input, $its1Output, $its2Output, $combinedOutput, $gbIndex) = @_;
###	my ($fastaSeqObj, $formattedSeqObj, $header, $desc);
###	my $in = Bio::SeqIO->new("-file"=>"<$input","-format"=>"Fasta");	#open all-seqs fasta file
###	my $its1Out = Bio::SeqIO->new("-file" => ">$its1Output", "-format" => "Fasta");
###	my $its2Out = Bio::SeqIO->new("-file" => ">$its2Output", "-format" => "Fasta");
###	my $combinedOut = Bio::SeqIO->new("-file" => ">$combinedOutput", "-format" => "Fasta");
###
###
###	my ($ITS1count, $ITS2count, $combinedCount) = (0,0,0);
###
###	while (defined($fastaSeqObj = $in->next_seq())) { #iterate over all the sequences
###
###		$header = getHeader($fastaSeqObj);
###		#my ($GI) = ($header =~ m/gi\|(\d+)\|/);
###		my ($GI) = ($header =~ m/gi\|([A-Z{2}\d]+\.\d+)\|/);
###		$desc = getFeaturesInfo($GI, $gbIndex);
###		$desc = formatKeyWords($desc);
###
###		#format the header to remove unnecessary ">" characters
###		$header =~ s/<(.+)>/[$1]/g;
###		$header =~ s/(>|<)//g;
###		my ($a, $b) = ($header =~ m/^(.+\|description\|)(.+)$/);
###
###		if($desc =~ m/ITS1/ && $desc =~ m/ITS2/){
###			#$header = "$a"."[1+2]"."$b"; #mark the sequence as a combined one
###			my ($id, $desc1) = ($header =~ /([^ ]*) (.*)/); #split the header into id and description (by the bioperl's genbank format): id = the header until the first space, desc = the rest of the header
###			$formattedSeqObj = Bio::Seq->new(-id => $id, -desc  => $desc1, -seq => $fastaSeqObj->seq() );
###
###			$combinedOut->write_seq($formattedSeqObj);
###			$combinedCount++;
###		}
###		elsif($desc =~ m/ITS1/ && $desc !~ m/ITS2/){
###			#$header = "$a"."[1]"."$b"; #mark the sequence as containing ITS1 only
###			my ($id, $desc1) = ($header =~ /([^ ]*) (.*)/); #split the header into id and description (by the bioperl's genbank format): id = the header until the first space, desc = the rest of the header
###			$formattedSeqObj = Bio::Seq->new(-id => $id, -desc  => $desc1, -seq => $fastaSeqObj->seq() );
###			$its1Out->write_seq($formattedSeqObj);
###
###			$ITS1count++;
###		}
###		elsif($desc !~ m/ITS1/ && $desc =~ m/ITS2/){
###			#$header = "$a"."[2]"."$b"; #mark the sequence as containing ITS2 only
###			my ($id, $desc1) = ($header =~ /([^ ]*) (.*)/); #split the header into id and description (by the bioperl's genbank format): id = the header until the first space, desc = the rest of the header
###			$formattedSeqObj = Bio::Seq->new(-id => $id, -desc  => $desc1, -seq => $fastaSeqObj->seq() );
###
###			$its2Out->write_seq($formattedSeqObj);
###			$ITS2count++;
###		}
###	}
###
###	return ($ITS1count, $ITS2count, $combinedCount);
###}


#def py_ITS_getFeaturesInfo():
#
#
#	return
#	sub getFeaturesInfo {
#	my ( $GI, $gbIndex ) = @_;
#	my ( @features, $featuresInfo, $featureValue, $featureType );
#	my $gbSeqObj = $gbIndex->fetch($GI);    # get the entry from GenBank
#	@features =	  $gbSeqObj->get_SeqFeatures;    #get the sequence's features (cds, source, gene....)
#
#	foreach my $featObj (@features) {
#		$featureValue = getFeatureValue($featObj);
#		$featureType  = $featObj->primary_tag;
#
#		#print "$featureType: $featureValue\n";
#		if ( !defined $featuresInfo ) {
#			$featuresInfo = $featureValue;
#		}
#		else {
#			$featuresInfo .= " | " . $featureValue;
#		}
#	}
#
#	return $featuresInfo;
#}

def merge_to_large_species(context):

	all_seqs = list(SeqIO.parse(context.fasta_seq_filename, "fasta"))

	for seq in all_seqs:
		tax_id = getPropertyFromFastaSeqHeader(seq.description, "taxonid")

		if tax_id in context.largeSpeciesTaxIDsList:
			logger.debug(seq.description)
			new_text = re.sub(r'(?s)(original_taxonid|)(.*?)(|seqid)', r"\1\3", seq.description)
			logger.debug(seq.description)
			return
#
		#parent_tax_id = context.parent_species_tax_id_by_subsp_tax_id[original_tax_id]
		#parent_name = context.parent_species_names_by_parent_tax_id[parent_tax_id]
		#f_subsp_merge_data.write(str(original_tax_id)+','+ str(original_name)+','+ str(parent_tax_id)+','+ str(parent_name)+'\n')

	return

#-------------------------------------------------------------------------------------------------------------
def create_ConstraintTree(working_dir,dict_tax_vs_Name,user_GRP_selection):

	ncbi = NCBITaxa()
	taxa_Grp_dict = {}

	# EXAMPLE for lineage:
	#>> ncbi.get_lineage(721805)
	#[1, 131567, 2759, 33090, 35493, 131221, 3193, 58023, 78536, 58024, 3398, 1437183, 71240, 91827, 1437201, 71275,
	# 91835, 3744, 3745, 171637, 721805]

	taxIDconstraint_dict = dict()
	dict_ID_Lineage_Grp=dict()
	dict_TaxID_lineage=dict()
	GroupTaxID_Const_dict=dict()
	#List of Constraint TaxID
	Constraint_TaxID_list=[]

	if user_GRP_selection == 'TaxaList':
		#read the TaxID for the constraint:
		with open(working_dir+'/ConstraintTaxIdList_user.txt', 'r') as const_taxIDList:
			taxa_cont=0
			for line in const_taxIDList:
				line_str = line.rstrip()
				if line_str != '':
					constr_taxId = int(line_str)
					logger.debug("constr_taxId: %d" %constr_taxId)
					taxIDconstraint_dict[constr_taxId] = taxa_cont
					if not Constraint_TaxID_list:
						Constraint_TaxID_list = [constr_taxId]
					else:
						Constraint_TaxID_list.append(constr_taxId)
					taxa_cont+=1
		#For each TaxID in the constraint list check parent and kid in list:
		# ncbi.get_lineage(int(ConstraintTax))
		# del checked_lin[-1]  # to remove ConstraintTax
		# list(reversed(checked_lin))
		RemovedFromSearch_list=[]
		spcParents_dict=dict()
		init_spc_list = list(dict_tax_vs_Name.keys())
		while init_spc_list:
			for spcTaxId in dict_tax_vs_Name.keys():
				if spcTaxId in RemovedFromSearch_list:
					continue
				else:
					RemovedFromSearch_list.append(spcTaxId)
					init_spc_list.remove(spcTaxId)
				spc_lineage = ncbi.get_lineage(int(spcTaxId))
				del spc_lineage[-1]  # to remove ConstraintTax
				revers_linage = list(reversed(spc_lineage))
				logger.debug("TaxID %d lineage: %s" %(int(spcTaxId), revers_linage))
				#Serach in species lineage:
				for lingTaxa in revers_linage:
					#logger.debug("Check if lingTaxa %d in Constraint list" % lingTaxa )
					if lingTaxa in Constraint_TaxID_list:
						spcParents_dict[spcTaxId]=lingTaxa

			for taxa in  dict_tax_vs_Name.keys():
				if taxa in spcParents_dict.keys():
					continue
				else:
					init_spc_list.append(taxa)
			logger.debug("spcParents_dict")
			logger.debug(spcParents_dict)
			logger.debug("RemovedFromSearch_list")
			logger.debug(RemovedFromSearch_list)
			logger.debug("dict_tax_vs_Name.keys()")
			logger.debug(dict_tax_vs_Name.keys())

		sys.exit()





		# Need to create a type of bubble sort for sorting taxIDs according to rank
		# For each one check linage and then take the lowest in its linsage check if in list and if so replace location if higher
		#for idConst in taxIDconstraint_dict.keys():
		#	id_lineage = ncbi.get_lineage(int(idConst))
		#	del id_lineage[-1]
		#	dict_TaxID_lineage[idConst] = id_lineage
#		#	logger.debug("Lineage of TaxID %d" %int(idConst))
		#	logger.debug("id_lineage)
		#	logger.debug(id_lineage)
#
		#for taxId_onTree in dict_tax_vs_Name:
		#	taxLineage = ncbi.get_lineage(int(taxId_onTree))
		#	del taxLineage[-1]
		#	taxLineage_rev = list(reversed(taxLineage))
		#	for lineageTaxa in taxLineage_rev:
		#		if lineageTaxa in Constraint_TaxID_list:
		#			taxa_Grp_dict[dict_tax_vs_Name[taxId_onTree]] = lineageTaxa	#Species belong to group#
#
		#logger.debug("taxa_Grp_dict")
		#logger.debug(taxa_Grp_dict)
#
		#list_groups_done=[]
		#for grp_Name in taxa_Grp_dict.values():
		#	if grp_Name in list_groups_done:
		#		print("Group already done %s" % grp_Name)
		#		continue
		#	else:
		#		list_groups_done.append(grp_Name)
		#	str_grp = '('
		#	for spc in taxa_Grp_dict.keys():
		#		str_grp+='%s,' %spc
		#	str_grp+=')'
		#	logger.debug(str_grp)
#
#
#
		##Sort TaxIds in order (rank wise):
		#for idConst in  dict_TaxID_lineage.keys():
		#	for idConstAll in dict_TaxID_lineage.keys():
		#		#check if in main TaxId lineage -> if so switch locations:
		#		if idConst in dict_TaxID_lineage[idConstAll]:
		#			idx_save = taxIDconstraint_dict[idConst] #save original taxId index location
		#			#Switch location:
		#			taxIDconstraint_dict[idConst] = taxIDconstraint_dict[idConstAll]
		#			taxIDconstraint_dict[idConstAll] = taxIDconstraint_dict[idConst]
#
#
		#logger.debug("Constraint TaxID after SORT:")
		#logger.debug(taxIDconstraint_dict)
#
		#tree = ncbi.get_topology([721813,3754,721805,721811,703251,180131,180133,4022])
		#tree.write(format=8, outfile=context.working_dir+'/ConstraintTree1.txt')
#
		#tree = ncbi.get_topology(Constraint_TaxID_list)
		#tree.write(format=8, outfile=context.working_dir + '/ConstraintTree2.txt')

		#

		#for idConst in  dict_TaxID_lineage.keys():
		#	for idConstAll in dict_TaxID_lineage.keys():
		#		if idConst == idConstAll:
		#			continue
		#		else:
		#			logger.debug("TaxID %s in TaxID %s" %(idConstAll,idConst))
		#			if idConst in dict_TaxID_lineage[idConstAll]:
		#				if idConst in dict_ID_Lineage_Grp:
		#					dict_ID_Lineage_Grp[idConst].append(idConstAll)
		#				else:
		#					dict_ID_Lineage_Grp[idConst]=[idConstAll]
		#logger.debug("dict_ID_Lineage_Grp")
		#logger.debug(dict_ID_Lineage_Grp)

		#For TaxID in keys there are groups that should be together.
		#TaxIds that are not in the keys/values are seperated groups:
		#logger.debug(dict_ID_Lineage_Grp)

		#Now need to go ove all the taxids in the species list and check which belongs to the groups of TaxIDs:
		#for spc in dict_tax_vs_Name.keys():
		#	if dict_ID_Lineage_Grp


		sys.exit()
	else:
		#Taxonomic constraint (wither Genus/Family/Both):
		for spc in dict_tax_vs_Name.keys():
			linage_list = ncbi.get_lineage(int(spc))
			# print(linage_list)
			linage_ranks = ncbi.get_rank(linage_list)
			for key in linage_ranks:
				if linage_ranks[key].lower()  == user_GRP_selection.lower():
					taxa_Grp_dict[spc] = str(key)

	# These are the grouping dictionaries for each taxon in the given list - either Genus or Family (user's choice)
	dict_name_grpName = dict()
	dict_numInGrp = dict()
	for key in taxa_Grp_dict:
		print(key)
		print("Grp %s -  is %s" % (taxa_Grp_dict, taxa_Grp_dict[key]))
	for value in taxa_Grp_dict.values():
		if value not in dict_numInGrp.keys():
			dict_numInGrp[value] = 1
		else:
			dict_numInGrp[value] += 1

	print(dict_numInGrp)

	with open(working_dir+'/ConstraintTree_user.txt','w') as const_f:
		str_constriant = ''
		str_constriant += ('(')
		list_groups_done = []
		for grp_Name in taxa_Grp_dict.values():
			number_of_species_in_grp = dict_numInGrp[grp_Name]
			if grp_Name in list_groups_done:
				print("Group already done %s" % grp_Name)
				continue
			else:
				list_groups_done.append(grp_Name)
			print(number_of_species_in_grp)
			if number_of_species_in_grp == 1:
				continue
			str_constriant += ('(')
			count = 1
			for key in taxa_Grp_dict.keys():
				name = dict_tax_vs_Name[key]
				if taxa_Grp_dict[key] == grp_Name:
					if count == number_of_species_in_grp:
						str_constriant += ('%s),' % name.replace(' ', '_'))
					else:
						str_constriant += ('%s,' % name.replace(' ', '_'))
						count += 1

		str_constriant += (');')
		str_constriant = str_constriant.replace(',);', ');')
		const_f.write(str_constriant)
	# Create dictionaries - for each group:

	return

#-------------------------------------------------------------------------------------------------------------
#This function will generate a constraint tree according to both Genera and Family:
def Create_Constraint_Both(working_dir,dict_tax_vs_Name):

	#For each family need to group all genera groups.
	user_GRP_selection = 'genus'
	ncbi = NCBITaxa()
	taxa_Grp_dict = {}
	for spc in dict_tax_vs_Name.keys():
		linage_list = ncbi.get_lineage(int(spc))
		linage_ranks = ncbi.get_rank(linage_list)
		for key in linage_ranks:
			if linage_ranks[key] == user_GRP_selection:
				taxa_Grp_dict[spc] = str(key)

	#These are the grouping dictionaries for each taxon in the given list - either Genus or Family (user's choice)
	dict_numInGrp = dict()
	for value in taxa_Grp_dict.values():
		if value not in dict_numInGrp.keys():
			dict_numInGrp[value]=1
		else:
			dict_numInGrp[value]+=1

	dict_GeneraID_vs_Groups = dict()
	list_groups_done=[]
	for grp_Name in taxa_Grp_dict.values():
		str_constriant = ''
		number_of_species_in_grp = dict_numInGrp[grp_Name]
		if grp_Name in list_groups_done:
			continue
		else:
			list_groups_done.append(grp_Name)
		if number_of_species_in_grp is '1':
			continue
		str_constriant += ('(')
		count = 1
		for key in taxa_Grp_dict.keys():
			name = dict_tax_vs_Name[key]
			if taxa_Grp_dict[key] == grp_Name:
				if count == number_of_species_in_grp:
					str_constriant += ('%s)' % name.replace(' ', '_'))
				else:
					str_constriant += ('%s,' % name.replace(' ', '_'))
					count += 1
		dict_GeneraID_vs_Groups[grp_Name] = str_constriant

	#Generate the family groups according to genera list:
	genera_list = list(dict_GeneraID_vs_Groups.keys())
	print(genera_list)
	user_GRP_selection = 'family'
	ncbi = NCBITaxa()
	taxa_FamilyGrp_dict = {}
	for spc in genera_list:
		linage_list = ncbi.get_lineage(int(spc))
		linage_ranks = ncbi.get_rank(linage_list)
		for key in linage_ranks:
			if linage_ranks[key] == user_GRP_selection:
				taxa_FamilyGrp_dict[spc] = str(key)
	dict_numInFamilyGrp = dict()
	for value in taxa_FamilyGrp_dict.values():
		if value not in dict_numInFamilyGrp.keys():
			dict_numInFamilyGrp[value]=1
		else:
			dict_numInFamilyGrp[value]+=1

	#Print the final constraint tree oncluding the Family grouping:
	dict_FamilyID_vs_Groups = dict()
	Final_GeneraFamily_str = '('
	with open(working_dir+'/ConstraintTree_user.txt','w') as const_f:
		list_groups_done = []
		for grp_Name in taxa_FamilyGrp_dict.values():
			str_constriant = '('
			number_of_elements_in_grp = dict_numInFamilyGrp[grp_Name]
			if grp_Name in list_groups_done:
				continue
			else:
				list_groups_done.append(grp_Name)
			print(number_of_elements_in_grp)
			if number_of_elements_in_grp is '1':
				continue
			count = 1
			for key in taxa_FamilyGrp_dict.keys():
				if taxa_FamilyGrp_dict[key] == grp_Name:
					if count == number_of_elements_in_grp:
						str_constriant += ('%s),' % dict_GeneraID_vs_Groups[key])
					else:
						str_constriant += ('%s,' % dict_GeneraID_vs_Groups[key])
						count += 1
			dict_FamilyID_vs_Groups[grp_Name] = str_constriant
			Final_GeneraFamily_str+=str_constriant

		Final_GeneraFamily_str+=');'
		Final_GeneraFamily_str = Final_GeneraFamily_str.replace(',);', ');')
		const_f.write(Final_GeneraFamily_str)
	print(dict_FamilyID_vs_Groups)

	return
#---------------------------------------------------------------------------------------------------------------
#In case the user choose to perfrom new clustering method in rerun Mode
def Rerun_Clustering(context):

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
		statusFail_LogFile(context, 'Less than 5 taxa were selected for each cluster.')
		return

	return

# -------------------------------------------------------------------------------------------------------
#This function will return the accession number of the sequence according to ID:
def get_clusterNum_from_path(path_str):

	#get accesion id from line like: gi|AF318735.1|taxonid|151425|organism|Prunus
	r = re.compile('concat\/(.*?)\/seqs-organism-concat')
	m = r.search(path_str)
	if m:
		return m.group(1)
	else:
		logger.debug("FAILED to find cluster number in path")
#-------------------------------------------------------------------------------------------------------