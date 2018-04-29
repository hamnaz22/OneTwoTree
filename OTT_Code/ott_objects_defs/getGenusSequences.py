import codecs
from Bio import SeqIO
import glob
import pandas as pd

from ott_objects_defs.ploidbCommon import *
from ott_objects_defs.ploidbUtils import turn_list_to_DB_list


__author__ = 'moshe'

MAX_SEQ_LENGTH_VITIS = 7000

#Get sequences from Database: Added Nov 2017 (adjust writeTaxonSequences_CdHit_KeepLargeData)
#
def DB_getGBsequences(context, fasta_filename, gb_filename, nextTaxId):

	genbank_db = ott_config['general']['DB_DIR'] + ott_config['concat_headir']['GENBANK_DB']
	max_seq_length = int(ott_config['general']['max_seq_length'])
	hugeTree_flag = ott_config['general']['Huge_trees']
	grp_list = ott_config['general']['GENBANK_GRP_LIST'].split(',')
	logger.debug("Using genabnk group list: %s" % grp_list)

	#Connect to Database to collect sequence data:
	conn = sqlite3.connect(genbank_db)
	curs = conn.cursor()
	list_for_query = turn_list_to_DB_list(context.species_list_ids)
	out_query = ''

	#Initial check for the number of sequences per TaxID -> this will save time later when collecting seqs for each taxID:
	#Create a dictionary for number of seqs per TaxID to enable the sepertion of large TaxIDs:
	#Check relevant groups and search all tables accordingly:
	count_grp = 1
	for grp_name in grp_list:
		grp_db_table = 'Genbank_seqs_' +grp_name
		out_query += "SELECT * from %s where TaxonID In (%s)" % (grp_db_table,list_for_query)
		if count_grp < len(grp_list):
			out_query += " UNION "
			count_grp+=1
	logger.debug("Perform %s" % (out_query))
	curs.execute(out_query)
	rows = curs.fetchall()

	#Make sure some data was collected: else notify the user:
	if not rows:
		logger.debug("No Data was fount in out Database, update your database and try again")
		with open(context.final_status, 'w') as f_status:
			f_status.write("No data found, update your genbank database")
		return


	#Insert data into csv file:
	headers_names_list = ["gb_group", "seqs_count", "TaxonId", "Accesion", "organism", "seq_length", "definition",
						  "gene_name_location", "voucher", "isolate", "product", "mol_type", "note", "desc",
						  "sequenceData"]
	df = pd.DataFrame(rows, columns=headers_names_list)
	f_csv_db = open(context.working_dir+'/GenbankDataFile.csv','w')
	df.to_csv(f_csv_db, header=False)
	#Close db connection:
	conn.close()
	f_csv_db.close()

	#Save data into context dictionaries/Lists .astype(str)
	context.taxon_Id_df = df['TaxonId'].astype(str)
	context.Accession_df = df['Accesion'].astype(str)
	context.organism_df = df['organism'].astype(str)
	context.seq_length_df = df['seq_length'].astype(str)
	context.definition_df = df['definition'].astype(str)
	context.gene_name_location_df = df['gene_name_location'].astype(str)
	context.desc_df = df['desc'].astype(str)
	context.mol_type_df = df['mol_type'].astype(str)
	context.note_df = df['note'].astype(str)
	context.product_df = df['product'].astype(str)
	context.sequenceData_df = df['sequenceData'].astype(str)

	logger.debug(len(context.taxon_Id_df))
	logger.debug(len(context.Accession_df))

	if len(context.taxon_Id_df) == 0:
		logger.debug("No TaxId was found in our Database")
		with open(context.final_status,'w') as f_stat:
			f_stat.write("Failed - The requested TaxIds are not included in our genabnk")

	df_count=0
	df_length=len(context.taxon_Id_df)
	cd_Dir = context.working_dir + '/CD_Hit_dir'
	create_dir_if_not_exists(cd_Dir)

	files_taxId_list =[]
	seq_count_perTaxID = dict()
	with open(cd_Dir + '/temp_AllSeqs.fasta', 'w') as t_all_temp:
		while df_count < df_length:
			taxon = str(context.taxon_Id_df[df_count])
			desc = str(context.desc_df[df_count])
			length = int(context.seq_length_df[df_count])
			seq = str(context.sequenceData_df[df_count])

			if taxon not in seq_count_perTaxID.keys():
				seq_count_perTaxID[taxon]=1
			else:
				seq_count_perTaxID[taxon]+=1

			f_tax = open(cd_Dir+'/'+'TaxID_'+str(taxon)+'_ForCdHit.fasta','a')
			files_taxId_list.append(f_tax)
			if length < max_seq_length:
				f_tax.write(desc + '\n' + seq + '\n')
				f_tax.close()

				no_of_seq_for_taxonid = seq_count_perTaxID[taxon]
				#Handle Large taxIds: more than 1000 seqs for one taxID
				if hugeTree_flag and no_of_seq_for_taxonid >= 1000:
					create_dir_if_not_exists(context.largeSpeciesDir)
					# Create ITS fasta file for Large Species ITS seqs that will be later added to the reg ITS file:
					ITS_large_fasta = open(context.largeSpeciesITS, 'w')

					f_large_log = open(context.largeSpeciesFileList, 'a')
					if taxon not in context.large_taxid_seqs_number_dict.keys():
						removed_f = open(context.summary_dir+'/RemovedTaxa.txt','a')
						removed_f.write("%s,Over 1000 seqs" %taxon)
						context.large_taxid_seqs_number_dict[taxon] = no_of_seq_for_taxonid

			df_count+=1

	#Move Large taxIds seq file to directory:

	for taxa in context.large_taxid_seqs_number_dict.keys():
		#taxa=str(taxa)
		large_fasta_filename = context.largeSpeciesDir + "/seqs_for_LargeTaxId_" + taxa + ".fasta"
		logger.debug("Copy seq file for Large taxa %s to %s" %(taxa,large_fasta_filename))
		shutil.copyfile(cd_Dir+'/'+'TaxID_'+str(taxa)+'_ForCdHit.fasta',large_fasta_filename)
		os.remove(cd_Dir+'/'+'TaxID_'+str(taxa)+'_ForCdHit.fasta')
		context.largeSpeciesTaxIDsList.append(taxa)
		f_large_log.write(large_fasta_filename + '\n')

		# remove ITS seqs from file:
		fasta_only_its_temp = large_fasta_filename + '_its_temp'
		remove_none_its_command = "perl " + context.cluster_script_dir + "/clustering/getSeqsWithITS.pl " + large_fasta_filename + " " + fasta_only_its_temp
		logger.info("Calling getSeqsWithITS (for ITS clustering) - " + remove_none_its_command)
		retval = exec_external_command_redirect_output(command_to_exec=remove_none_its_command,
													   outfile=context.largeSpeciesDir + "/getSeqsWithITS.out",
													   errfile=context.largeSpeciesDir + "/getSeqsWithITS.err")
		# add to main large ITS file:
		with open(fasta_only_its_temp, 'rU') as itsTemp_file:
			seq_records = SeqIO.parse(itsTemp_file, 'fasta')
			SeqIO.write(seq_records, ITS_large_fasta, 'fasta')

	if os.path.exists(context.largeSpeciesFileList):
		f_large_log.close()

	#Once we have all files we need to go over all of them CD-hit them and combine to one file:
	no_of_seq_for_taxonid=0
	files = glob.glob(cd_Dir+'/*_ForCdHit.fasta')
	logger.debug(len(files))
	for Seq_file in files:
		shutil.copyfile(Seq_file,context.fasta_cdhit_temp_input_f)
		#logger.info('Running CD-HIT on Taxid - %s' % nextTaxId)
		with open(context.fasta_cdhit_temp_input_f, 'r') as cd_input_f:
			lines_input = cd_input_f.readlines()
			num_of_seq_BeforeCdHit = (len(lines_input) / 2)
		#Check if number of seqs is less than 5 -> copy file as is:
		if num_of_seq_BeforeCdHit <= 5:
			shutil.copyfile(context.fasta_cdhit_temp_input_f,context.fasta_cdhit_temp_output_f)
		else:
			CDHITcommand_to_exec = 'cd-hit -i %s -d 0 -o %s -c 0.98 -n 5 -G 1 -g 1 -b 25' % (context.fasta_cdhit_temp_input_f,context.fasta_cdhit_temp_output_f)
			exec_external_command_redirect_output(command_to_exec=CDHITcommand_to_exec,outfile=context.working_dir + "/CD-HIT_perTaxID.out",errfile=context.working_dir + "/CD-HIT_perTaxID.err")
		if fasta_filename is not None:
			with open(fasta_filename, 'a') as f1:
				read_cd_File = open(context.fasta_cdhit_temp_output_f, 'r')
				lines = read_cd_File.readlines()
				for line in lines:
					if line.strip(): f1.write(line) # write only non-blank lines
			num_of_seq_AfterCdHit = (len(lines)/2)
			no_of_seq_for_taxonid+=num_of_seq_AfterCdHit
			logger.debug("Wrote " + str(num_of_seq_AfterCdHit) + ' out of ' + str(num_of_seq_BeforeCdHit) + " new sequences in fasta format (after CD-Hit)")



	temp_fasta_file = open(context.fasta_cdhit_temp_input_f, 'w')
	temp_fasta_o_file = open(context.fasta_cdhit_temp_output_f, 'w')
	temp_fasta_file.close()
	temp_fasta_o_file.close()

	# Check if hugeTree flag:
	hugeTree_flag = ott_config['general']['Huge_trees']


	#sequences_path = getSpeciesSequenceFileLocation(nextTaxId)  # getting gene data from the files downloaded from genebank
	return no_of_seq_for_taxonid


def writeTaxonSequences(fasta_filename, gb_filename, nextTaxId):
	sequences_path = getSpeciesSequenceFileLocation(nextTaxId)  # getting gene data from the files downloaded from genebank
	if sequences_path is None:
		# TODO: add code that checks if this is legitimate (i.e by querying NCBI online) or not
		# This can happen in cases where the taxon itself doesn't have sequences but its subspecies do
		logger.warning("Path not found for taxon " + str(nextTaxId) + ". Skipping this taxon")
		return 0
	else:
		logger.info("getting sequences for " + nextTaxId + " from " + str(sequences_path))
		speciesSeqRecords = SeqIO.parse(sequences_path, "genbank")
		# The data is sorted so that on multiple runs of the pipeline, it would be easier to use cached results (which are based on file CRC)
		speciesSeqList = get_ordered_seqs(list(speciesSeqRecords))
		no_of_seq_for_taxonid = len(speciesSeqList)
		logger.debug("found " + str(no_of_seq_for_taxonid) + " sequences for " + nextTaxId)

		speciesSeqListWrittenToFasta = list()
		if fasta_filename is not None:
			all_genus_seq_fasta_handle = open(fasta_filename, "a")

			mac_seq_length = int(ott_config['general']['max_seq_length'])

			for seq in speciesSeqList:
				if len(seq.seq) > mac_seq_length:
					logger.info("Seq length is %d for seq %s - ignoring this seq" % (len(seq.seq), seq.description))
				# Since Vitis
				# elif genusName.lower() == 'vitis' and (
				# len(seq.seq) > MAX_SEQ_LENGTH_VITIS or "whole genome shotgun" in seq.description):
				# logger.info("Seq length is %d for seq %s - ignoring this seq" % (len(seq.seq), seq.description))
				else:
					writeSequenceInFastaFormat(all_genus_seq_fasta_handle, nextTaxId, seq)
					speciesSeqListWrittenToFasta.append(seq)

			all_genus_seq_fasta_handle.close()
		# logger.debug("Wrote " + str(no_of_seq_for_taxonid) + " new sequences in fasta format")

		if gb_filename is not None:
			# logger.debug("Writing to GB file")
			all_genus_seq_gb_handle = open(gb_filename, "a")
			# TODO: check the sequece length and filter if its too long
			SeqIO.write(speciesSeqListWrittenToFasta, all_genus_seq_gb_handle, "genbank")
			all_genus_seq_gb_handle.close()
		# logger.debug("Wrote " + str(no_of_seq_for_taxonid) + " new sequences in genbank format")

		return no_of_seq_for_taxonid


def Create_ListofSpecies(context):
	species_list=list()
	for name in context.taxa_list:
		name_clean=name.replace('\'', '')
		name_clean=name_clean.replace('\"', '')
		name_clean=name_clean.strip()
		print(name_clean)
		new_name = "'"+name_clean+"'"
		print(name)
		species_list.append(new_name)
	species_list_str = ','.join(species_list)
	return species_list_str


# getting taxa sequences from taxa id's and saving it to the fasta and gb files
def get_taxa_sequences(context,species_ids, fasta_filename, gb_filename, p):
	if os.path.exists(fasta_filename):
		os.remove(fasta_filename)

	if context.RerunID == 'None':
		if os.path.exists(gb_filename):
			os.remove(gb_filename)

	numberOfSequenceForSpecies = 0
	species_ids = sorted(species_ids)
	taxon_id='0'

	# adding seq of all taxons in list to .gb and .fasta files from the genebank files
	DB_StandAlone=1
	if DB_StandAlone == 0:
		for taxon_id in species_ids:
			# With or Without CD-Hit:
			CDHIT_FLAG = ott_config['general']['cd_hit_flag']
			print(CDHIT_FLAG)
			if CDHIT_FLAG == 'FALSE':
				numberOfSequenceForSpecies += writeTaxonSequences(fasta_filename, gb_filename, taxon_id)
			else:
				numberOfSequenceForSpecies += DB_getGBsequences(context, fasta_filename, gb_filename, taxon_id)
				#numberOfSequenceForSpecies += writeTaxonSequences_CdHit_KeepLargeData(context, fasta_filename, gb_filename, taxon_id)
	else:
		numberOfSequenceForSpecies += DB_getGBsequences(context, fasta_filename, gb_filename, taxon_id)

	logger.debug("Saved to " + fasta_filename)
	#remoce cd-hit files:
	if os.path.exists(context.working_dir + "/CD-HIT_perTaxID.out"): os.remove(context.working_dir + "/CD-HIT_perTaxID.out")
	if os.path.exists(context.working_dir + "/CD-HIT_perTaxID.err"): os.remove(context.working_dir + "/CD-HIT_perTaxID.err")

	return numberOfSequenceForSpecies
