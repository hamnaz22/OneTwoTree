#--------------------------------------------------------------------------------------------
#--  How to use me:
#--	 Make sure you run it from the main OTT directory (where you extracted the OTT package)
#--  Running cmd: python GenBank_OTT.py pln-pri-mam
#--------------------------------------------------------------------------------------------
import sys
import os
import csv
import pandas as pd
import numpy as np
import sqlite3
import glob
from Bio import SeqIO


Debug_Flag = 1

#-----------------------------     DOWNLOAD smp NCBI files (nodes and names)    -----------------------------------
def Create_download_file(path_for_download):

	dmp_Dir = path_for_download+'ncbi_dmp/'
	if not os.path.exists(dmp_Dir):
		os.makedirs(dmp_Dir)
	os.system('cd %s' %dmp_Dir)
	with open(dmp_Dir+'down_dmp.sh','w') as f_down_dmp:
		f_down_dmp.write("#!/bin/tcsh\n\n")			# #!/bin/bash
		f_down_dmp.write("#$ -N down_dmp\n")
		f_down_dmp.write("#$ -S /bin/tcsh\n")
		f_down_dmp.write("#$ -cwd\n")
		f_down_dmp.write("\n")
		f_down_dmp.write("module load python/anaconda3-5.0.0\n\n")
		f_down_dmp.write("cd %s\n" %dmp_Dir)
		f_down_dmp.write("ftp -n -i ftp.ncbi.nih.gov<<EOF\n")
		f_down_dmp.write("user anonymous michal.md75@gmail.com\n")
		f_down_dmp.write("cd pub/taxonomy/\n")
		f_down_dmp.write("mget taxdump.tar.gz\n")
		f_down_dmp.write("quit\n")
		f_down_dmp.write("EOF\n")
		f_down_dmp.write("\n")
#
		f_down_dmp.write("cd %s\n" %dmp_Dir)
		f_down_dmp.write("tar -xvzf taxdump.tar.gz\n\n")
#
	f_down_dmp.close()

	os.system('tcsh %s' % (dmp_Dir+'down_dmp.sh'))
	return 'done'

#-----------------------------     DOWNLOAD smp NCBI files (nodes and names)    -----------------------------------
def grp_genabnk_download(grp_name,path_for_download,f_log):

	grp_Dir = path_for_download+'gb_'+grp_name
	if not os.path.exists(grp_Dir):
		os.makedirs(grp_Dir)
		print("create dir %s" %grp_Dir)
	else:
		f_log.write("group %s dir already exist\n" % grp_name)
		return
	os.system('cd %s' %grp_Dir)
	with open(grp_Dir+'/'+grp_name+'_download.sh','w') as f_gb_download:
		f_gb_download.write("#!/bin/tcsh\n\n")			# #!/bin/bash
		f_gb_download.write("#$ -N GenBank-download\n")
		f_gb_download.write("#$ -S /bin/tcsh\n")
		f_gb_download.write("#$ -cwd\n")
		f_gb_download.write("\n")
		f_gb_download.write("cd %s\n" %grp_Dir)
		f_gb_download.write("ftp -n -i ftp.ncbi.nih.gov<<EOF\n")
		f_gb_download.write("user anonymous michal.md75@gmail.com\n")
		f_gb_download.write("cd genbank\n")
		f_gb_download.write("mget gb%s*\n" %grp_name)
		f_gb_download.write("quit\n")
		f_gb_download.write("EOF\n")
		f_gb_download.write("\n")

		f_gb_download.write("cd %s\n" %grp_Dir)
		f_gb_download.write("gunzip *\n\n")

	f_gb_download.close()

	os.system('tcsh %s' % (grp_Dir+'/'+grp_name+'_download.sh'))
	f_log.write("group %s download from genbank was completed\n" %grp_name)
	print("group %s download from genbank was completed" %grp_name)
	return 'done'


#---------------------------------    Create Nodes dictionaries    --------------------------------------
def Create_ParentsChildDict_db(ncbi_nodes_file):

	#Nodes_Parents_dict=dict()
	#Nodes_Rank_dict=dict()
	#Read files into dataFrame and conver to Dictionaries:
	df_1 = pd.read_csv(ncbi_nodes_file,sep='|',names =  ["TaxId", "ParentTaxId", "Rank", "A", "B", "C", "D", "E", "F",
														 "G", "H", "I", "J"],index_col=False).astype(str)
	df_2 = pd.read_csv(ncbi_nodes_file,sep='|',names =  ["TaxId", "ParentTaxId", "Rank", "A", "B", "C", "D", "E", "F",
														 "G", "H", "I", "J"],index_col=False).astype(str)

	Nodes_Parents_dict = df_1.set_index('TaxId')['ParentTaxId'].to_dict()
	Nodes_Rank_dict = df_2.set_index('TaxId')['Rank'].to_dict()

	return Nodes_Parents_dict,Nodes_Rank_dict


#-------------- Fast import data into DataBase   ------------------------------------------
def cvd_file_into_DB(arr_to_insert,db_file):

	conn = sqlite3.connect(db_file)
	curs = conn.cursor()
	curs.execute('''CREATE TABLE ncbi_names_db (Rank TEXT, TaxID INT, ParentTaxID INT, Name TEXT, Other TEXT, 
	Type TEXT)''')
	results_header = "Rank, TaxID, ParentTaxID, Name, Other, Type"
	curs.executemany("INSERT INTO ncbi_names_db (" + results_header + ") VALUES (?,?,?,?,?,?);", arr_to_insert)

	# Save (commit) the changes
	conn.commit()
	conn.close()
	return

#--------------------  Prepare data for DataBase
def prepare_for_DB(ncbi_names_file,nodes_Parents_dict,nodes_Rank_dict):

	ar_ar = []
	with open(ncbi_names_file, 'r') as ncbi_data_f:
		idx_cnt=0
		for line in ncbi_data_f:
			idx_cnt+=1
			if idx_cnt%100000 == 0:
				print(idx_cnt)
			parent_taxId = '-'
			line_parse = line.replace('\t', '')
			line_parse = line_parse.replace('\'', '')
			line_parse = line_parse.replace('\"', '')
			line_parse = line_parse.split('|')
			#Check if Type is relevant:
			if line_parse[3] == 'type material' or line_parse[3] == 'includes':
				continue
			if line_parse[0] in nodes_Parents_dict.keys():
				parent_taxId = str(nodes_Parents_dict[line_parse[0]])
				rank = str(nodes_Rank_dict[line_parse[0]]).replace('\t', '')
			ar_ar.append([rank, line_parse[0], parent_taxId, line_parse[1], line_parse[2], line_parse[3]])

	return ar_ar

#--------------------       ????????????????????????        ------------------------------------------
def return_seq_desc(seq,f_log):

	if Debug_Flag == 1: print("-------------------   %s  ----------------------" %seq.annotations["organism"])
	organism = seq.annotations["organism"]
	taxon_id = 'xxxxxx'
	isolate_name = '-'
	spec_voucher = '-'
	product_name = '-'
	mol_type = '-'
	note_name = '-'

	gene_names_location_dict=dict()
	for feature in seq.features:
		for key in feature.qualifiers.keys():
			if feature.type == 'gene':
				start = int(feature.location.start)
				stop = int(feature.location.end)
				str_start_stop=str(start) +'-'+str(stop)
				if "gene" not in feature.qualifiers.keys():
					gene_name ='NotClear'
				else:
					gene_name = feature.qualifiers["gene"][0]
				gene_names_location_dict[gene_name]=str_start_stop
			if key == 'db_xref':
				#To skip all other db_xref values, such as BOLD:
				for db_data in feature.qualifiers["db_xref"]:
					if 'taxon' in db_data:
						taxon_id = db_data.strip('taxon:')
			if key == 'isolate':
				isolate_name = feature.qualifiers["isolate"][0]
			if key == 'specimen_voucher':
				spec_voucher = feature.qualifiers["specimen_voucher"][0]
			if key == 'product':
				product_name = feature.qualifiers["product"][0]
			if key == 'mol_type':
				mol_type = feature.qualifiers["mol_type"][0]
			if key == 'note':
				note_name = feature.qualifiers["note"][0]

	gene_concat=''
	if gene_names_location_dict:
		for key in gene_names_location_dict.keys():
			gene_concat+=key +':'+gene_names_location_dict[key]+'|'
	else:
		gene_concat = '-'
	#isol_prod_mol_note = isolate_name + '|' + product_name+ '|' + mol_type+ '|' + note_name
	#gi = seq.annotations["accessions"][0] #accessions
	desciption_str = seq.description.replace("'","").replace('"','').replace(',','')
	organism = organism.replace("'","").replace('"','').replace(',',' ')
	seq_desc = ">gi|%s|taxonid|%s|organism|%s|seqid|%s|description|%s" %(seq.id, taxon_id,organism,seq.id,
																		 desciption_str)
	return seq_desc,taxon_id,seq.id,organism,desciption_str,gene_concat,spec_voucher,isolate_name,product_name,mol_type,\
		   note_name
#--------------------  Create Group table and insert data ------------------------------------------
def insert_grp_tableData(Seq_file_path,db_grp,grp_name,f_log):

	os.chdir(Seq_file_path)
	conn = sqlite3.connect(db_grp)
	curs = conn.cursor()
	#Create Table for group:
	grp_table_name = 'Genbank_seqs'+'_'+grp_name
	curs.execute(
		'''CREATE TABLE %s (gb_group STRING, seqs_count STRING, TaxonId STRING, Accesion STRING, organism STRING, 
		seq_length STRING, definition STRING, gene_name_location STRING, voucher STRING, isolate STRING, product STRING, 
		mol_type STRING, note STRING, desc STRING, sequenceData STRING, PRIMARY KEY(seqs_count,TaxonId))'''
		%grp_table_name)
	conn.commit()
	db_arr=[]
	seq_count_perTaxID=dict()
	for Seq_file in glob.glob("*.seq"):
		f_log.write("Read sequences from: %s\n" %Seq_file)
		print("Read sequences from: %s\n" %Seq_file)
		SeqRecords = SeqIO.parse(Seq_file, "genbank")
		SeqRecords_list = list(SeqRecords)
		SeqList = sorted(SeqRecords_list, key=lambda seq: seq.description)
		no_of_seq_for_taxonid = len(SeqList)
		max_seq_length = 10000
		spec_count=0
		for seq in SeqList:
			seq_desc,taxon_id,Accesion,Organism,definition,gene_concat,spec_voucher,isolate,product,mol_type,note = \
				return_seq_desc(seq,f_log)
			if taxon_id not in seq_count_perTaxID.keys():
				seq_count_perTaxID[taxon_id]=1
			else:
				seq_count_perTaxID[taxon_id]+=1
			seq_length = len(seq.seq)
			if seq_length <= max_seq_length:
				seq_str = seq.seq
				db_arr.append([grp_name,seq_count_perTaxID[taxon_id],taxon_id,Accesion,Organism,seq_length,definition,
							   gene_concat,spec_voucher,isolate,product,mol_type,note,seq_desc,str(seq_str)])
			spec_count+=1

	f_log.write("Insert data to Database\n")
	results_header = "gb_group,seqs_count,TaxonId,Accesion,organism,seq_length,definition,gene_name_location,voucher," \
					 "isolate,product,mol_type,note,desc,sequenceData"
	curs.executemany("INSERT INTO " +grp_table_name+" (" + results_header + ") VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);",
					 db_arr)
	conn.commit()
	conn.close()
	return

#--------------------------------------------------------------------------------------------
#--                                                                                        --
#------------------------------------   MAIN   ----------------------------------------------
#--                                                                                        --
#--------------------------------------------------------------------------------------------
grp_names = sys.argv[1]
path_for_download = sys.argv[2] #ploidb_config['general']['DB_DIR']
grp_list = []
if len(sys.argv) < 1:
	sys.exit("Please insert the group names you would like to download from NCBI (e.g. pln or pln|mam etc.")
for grp in grp_names.split('-'):
	grp_list.append(grp)
print (grp_list)

if os.path.exists(path_for_download):
	log_debug = path_for_download + 'LOG_debug.txt'
	ncbi_names_file = path_for_download + 'ncbi_dmp/names.dmp'
	ncbi_nodes_file = path_for_download + 'ncbi_dmp/nodes.dmp'
	NCBI_db_file = path_for_download + 'NCBI_names_NR.db'
	NCBI_dbFast_file = path_for_download + 'NCBI_NR.db'
	GenbankDB = path_for_download + 'GenbankDB_grp_test.db'

	#Download New dmp files:
	if Create_download_file(path_for_download) == 'done': #dmp files are ready
		Nodes_Parents_dict,Nodes_Rank_dict = Create_ParentsChildDict_db(ncbi_nodes_file)
		print("End of stage1")
		arr_to_insert = prepare_for_DB(ncbi_names_file,Nodes_Parents_dict,Nodes_Rank_dict)
		print("End of stage2")
		cvd_file_into_DB(arr_to_insert,NCBI_dbFast_file)
		print("The end: NCBI_NR database is ready :)")


	# Download Genabnk groups to Genabnk DB:
	f_log = open(log_debug,'w')
	for grp_name in grp_list:
		grp_genabnk_download(grp_name,path_for_download,f_log)
		print('Group %s downloaded successfully' %grp_name )
		insert_grp_tableData(path_for_download+'gb_'+grp_name,GenbankDB,grp_name,f_log)
		print('Insert group %s into DB %s' %(grp_name,GenbankDB))
		#remove group files:
		os.system("rm -rf %s" %(path_for_download+'gb_'+grp_name))

	f_log.close()

else:

	print("DB-dir : Database directory is missing...")
#--------------------------------------------------------------------------------------------
#--                                                                                        --
# 					          				THE END   	 									#
#--                                                                                        --
#--------------------------------------------------------------------------------------------
