import time
import pandas as pd
import os
import re
import glob
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from ott_objects_defs.ObjectJSONEncoder import ObjectJSONEncoder
from ott_objects_defs.PloiDbContext import PloiDbContext
from ott_objects_defs.ploidbCommon import *
from ott_objects_defs.handleMultAccessions import get_taxon_gis_key

from buildTaxaTree import perform_filter_msa
from pandas import Series


__author__ = 'ItayM3'

#-------------------------------------------------------------------------------------------------------
def prepare_files_for_concat(indexIts12FileName):

	with open(indexIts12FileName,'r') as f_list:
		for line in f_list:
			file_name_path = line.strip()
			in_fasta_filename = file_name_path
			out_organism_desc_fasta_filename = file_name_path+'_edit'
			rewrite_fasta_with_organism_name_as_desc(in_fasta_filename, out_organism_desc_fasta_filename, add_gi=False)
			shutil.copyfile(out_organism_desc_fasta_filename,file_name_path)
			logger.debug("File %s headers edited to be organism only" %file_name_path)

	return


#-------------------------------------------------------------------------------------------------------
def is_empty_file(f_handle):

	for line in f_handle:
		if not line.strip():
			continue
		else:
			#Not empty
			return 1
	#Empty File
	return 0
#-------------------------------------------------------------------------------------------------------
#                   ITS FLOW
#-------------------------------------------------------------------------------------------------------
# 1. pickOneITSTypePerSpeciesFasta_py: 	split sequence fasta file by species and locate each file in taxId directory
# 										under 'species_all' directory. for each file perform split by ITS type using
# 										- splitITS_py
# 		1.1 splitITS_py:
# 					create 3 ITS fasta files accoring to type -> combined, ITS1_only & ITS2_only.
# 					return the number of each type for this species.
# 					* for later use - calc_longest_accordingToType: calculate the length of the longest sequence from
# 						each type for later filtering process
# 		1.2 Choose representative:
# 					for each species choose a representative sequence:
#					If combined exist: use pickFromFile_py to select the best seq using blastall and best score
#					compared with other seqs. Copy this seq to file : combined.fasta
#					ELSE:
#					If ITS1 or ITS2: use pickFromFile_py to select and copy seqs to ITS1_only.fasta
# 					and ITS2_only.fasta
#
# 2. ITS_MSA_py / ITS_CLUSTALO_py:
#		2.1 If combined and SEP perform Adjust direction on both, else just on combined. Using adjustDirFasta.py on:
#			SEP_ITS1+ITS2.fasta: file created by concatenation of both ITS1_only.fasta & ITS2_only.fasta
#		2.2 Perform mafft (mafft --auto --ep 0) on combined.fasta -> Output to combined.msa
#		2.3 Add fragments: mafft --addfragments ../SEP_ITS1+ITS2.fasta --multipair ../combined.msa > ../combined+sep.msa
#		2.4 If NO combined:
#			Adjust direction amd mafft on ITS1
#			Adjust direction amd mafft on ITS2
#			Mafft each file and concatenate them using 'cat'.
#			Perform another mafft on concatenated file: mafft --auto --quiet %s/SEP_ITS1+ITS2_temp.msa ????? Check if needed?!!!!!
#
# 3. After MSA - selecting a single ITS seq per species	 - pickOneSeqPerSpeciesMSA:
#		3.1 take the Output msa file from stage 2 and split by species under 'species' directory again using
#			SplitFastaSeqsBySpecies.pl
#		3.2 If combined exists: add the aligned seq to Final_records
#		3.3 ELSE: if (ITS1count > 0 and ITS2count > 0): Merging(Combined Exists)/Appending ITS1 and ITS2:
#			Use getChosenSequence to perform append or merge on ITS1 and ITS2.
# 		3.4 ELSE if just one type add this seq as is (msa)
#		Final output is an aligned file: oneSeqPerSpecies.msa
#
# 4. Guidance:
# 		4.1 Prepare fasta file of the combined and merged seqs (this will be the bsae) and the rest of ITS1&2 that were
#			not merged will be added as fragments.
#		4.2 Run Guidacne mafft with addfragments
#		4.3 Remove the sequences marked by guidance from 'oneSeqPerSpecies.msa' and rerun mafft with add frag.
#
#
#-------------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------------------
# This function will run guidance ...
#--------------------------------------------------------------------------------------------------------------
def new_guidance_for_its(context,sequencesOfCoreMSA, fragments, outDir):

	if context.UserFlags_dict['MSA_Software'] == 'ClustalOmega':
		guidanceMSA_name = 'CLUSTALW'
	else:
		guidanceMSA_name = 'MAFFT'
	guidance_path = ott_config['diff_soft']['Guidance']

	#Check if needs to add fragments of spe1_2:
	if fragments != 'None':
		guidance_cmd = "perl %s --program GUIDANCE --seqFile %s --msaProgram MAFFT --seqType nuc --outDir %s " \
					   "--MSA_Param \"\\-\\-adjustdirection \\-\\-addfragments %s \\-\\-multipair\"" %(guidance_path,sequencesOfCoreMSA,outDir,fragments)
	else:
		guidance_cmd = "perl %s --program GUIDANCE --seqFile %s --msaProgram MAFFT --seqType nuc --outDir %s " \
					   "--MSA_Param \"\\-\\-adjustdirection\"" % (
					   guidance_path, sequencesOfCoreMSA, outDir)

	#This is the final msa that is copied into the ITS cluster dir under the concat dir:
	fasta_final_output = outDir + '/oneSeqPerSpecies.msa'
	logger.debug("ITS guidance cmd: %s" % guidance_cmd)
	retval = exec_external_command_redirect_output(guidance_cmd)

	#Once guidance ran we check if any seqs need to be removed:
	remove_seq_f =  outDir + "/Seqs.Orig.fas.FIXED.Removed_Seq.With_Names"
	removed_acc_id_list=[]
	if os.stat(remove_seq_f).st_size != 0:
		records_list = list(SeqIO.parse(remove_seq_f, "fasta"))
		records_core_list = list(SeqIO.parse(sequencesOfCoreMSA, "fasta"))
		records_frag_list = list(SeqIO.parse(fragments, "fasta"))
		for seq in records_list:
			seq_id = getPropertyFromFastaSeqHeader(seq.description, "gi")
			removed_acc_id_list.append(seq_id)
		for seq in records_core_list:
			seq_id = getPropertyFromFastaSeqHeader(seq.description, "gi")
			if seq_id not in removed_acc_id_list:
				Final_records.append(seq_id)
		SeqIO.write(Final_records, sequencesOfCoreMSA, "fasta")
		for seq in records_frag_list:
			seq_id = getPropertyFromFastaSeqHeader(seq.description, "gi")
			if seq_id not in removed_acc_id_list:
				Final_records.append(seq_id)
		SeqIO.write(Final_records, records_frag_list, "fasta")



	if retval != 0:
		with open(context.final_status, 'w') as f_status:
			f_status.write("ITS-guidance-Failed")
		raise Exception("ITS guidance failed: Failed to choose the final tree")
	else:
		last_output_from_guidance = outDir + "/MSA."+guidanceMSA_name+".Without_low_SP_Col.With_Names"
		logger.info("Copying final guidance/alignment results from %s to %s" % (last_output_from_guidance, fasta_final_output))
		shutil.copy(last_output_from_guidance, fasta_final_output)


	return

#--------------------------------------------------------------------------------------------------------------
# This function will run guidance either on combined+sep / combined files.
# for combined+sep it will use addfrag option. also the msa will be performed according to user choice
# Guidance output:
# MSA without problematic columns: MSA.MAFFT.Without_low_SP_Col.With_Names
# if 'Seqs.Orig.fas.FIXED.Removed_Seq' file is not empty:
# take sequence file(un-aligned): Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names
# run guidance again and then take: MSA.MAFFT.Without_low_SP_Col.With_Names
#--------------------------------------------------------------------------------------------------------------
def guidance_for_its(context,sequencesOfCoreMSA, fragments, outDir):

	#if context.UserFlags_dict['MSA_Software'] == 'ClustalOmega':
	#	guidanceMSA_name = 'CLUSTALW'
	#else:
	guidanceMSA_name = 'MAFFT'
	guidance_path = ott_config['diff_soft']['Guidance']

	if fragments != 'None':
		guidance_cmd = "perl %s --seqFile %s --program GUIDANCE --msaProgram MAFFT --seqType nuc --outDir %s " \
					   "--MSA_Param \"\\-\\-adjustdirection \\-\\-addfragments %s \\-\\-multipair\"" %(guidance_path,sequencesOfCoreMSA,outDir,fragments)
		fasta_final_output = outDir + '/combined+sep.msa'
	else:
		guidance_cmd = "perl %s --seqFile %s --program GUIDANCE --msaProgram MAFFT --seqType nuc --outDir %s " \
					   "--MSA_Param \"\\-\\-adjustdirection\"" % (
					   guidance_path, sequencesOfCoreMSA, outDir)
		fasta_final_output = outDir + '/combined.msa'


	logger.debug("ITS guidance cmd: %s" % guidance_cmd)
	retval = exec_external_command_redirect_output(guidance_cmd)#

	if retval != 0:
		with open(context.final_status, 'w') as f_status:
			f_status.write("ITS-guidance-Failed")
		raise Exception("ITS guidance failed: Failed to choose the final tree")
	else:
		last_output_from_guidance = outDir + "/MSA."+guidanceMSA_name+".Without_low_SP_Col.With_Names"
		logger.info("Copying final guidance/alignment results from %s to %s" % (last_output_from_guidance, fasta_final_output))
		shutil.copy(last_output_from_guidance, fasta_final_output)


	return


#-------------------------------------------------------------------------------------------------------
#This function will calc the longest seq for each type: its1/its2/combined:
def calc_longest_accordingToType(context):
	# calc avg length of each type:
	com_avg_length = 0;	comb_total = 0;	comb_count = 0
	its1_avg_len = 0; its1_total = 0; its1_count = 0
	its2_avg_len = 0; its2_total = 0; its2_count = 0
	for id in context.its_accession_ids:
		if id in context.its_accession_vs_type:
			if context.its_accession_vs_type[id] == 'combined':
				comb_total += context.its_accession_vs_length[id]
				comb_count += 1
			if context.its_accession_vs_type[id] == 'its1':
				its1_total += context.its_accession_vs_length[id]
				its1_count += 1
			if context.its_accession_vs_type[id] == 'its2':
				its2_total += context.its_accession_vs_length[id]
				its2_count += 1
	if comb_count != 0:
		context.its_type_avg_len['combined'] = comb_total / comb_count
	else:
		context.its_type_avg_len['combined'] = 0
	if its1_count != 0:
		context.its_type_avg_len['its1'] = its1_total / its1_count
	else:
		context.its_type_avg_len['its1'] = 0
	if its2_count != 0:
		context.its_type_avg_len['its2'] = its2_total / its2_count
	else:
		context.its_type_avg_len['its2'] = 0
	logger.debug("Summary of ITS type:\n")
	logger.debug("Combined count: %d, avg length: %f\n" % (comb_count, context.its_type_avg_len['combined']))
	logger.debug("ITS1 count: %d, avg length: %f\n" % (its1_count, context.its_type_avg_len['its1']))
	logger.debug("ITS2 count: %d, avg length: %f\n" % (its2_count, context.its_type_avg_len['its2']))

	return

#-------------------------------------------------------------------------------------------------------
#This function will return the accession number of the sequence according to ID:
def get_acc_from_id(acccession_id_line):

	#get accesion id from line like: gi|AF318735.1|taxonid|151425|organism|Prunus
	r = re.compile('gi\|(.*?)\|taxonid')
	m = r.search(acccession_id_line)
	if m:
		return m.group(1)
	else:
		logger.debug("FAILED to find accesion id in line (get_acc_from_id)")

# -------------------------------------------------------------------------------------------------------
#This function will return the id from accession_id_list of the seq with the most similar length of it's type:
def return_longest_seq(context,accession_id_list):

	min_distance = 999999
	selected_id = accession_id_list[0]
	for id in accession_id_list:
		id_type = context.its_accession_vs_type[id]
		id_length = context.its_accession_vs_length[id]
		if abs(float(id_length - context.its_type_avg_len[id_type])) < min_distance:
			logger.debug('return_longest_seq: selected length %d, Avg length %f' %(id_length,context.its_type_avg_len[id_type]))
			selected_id = id

	return selected_id



#create sequence object for the merged sequence created from seq1Obj and seq2Obj
def createMergedSeqObj_py(seq1Obj, seq2Obj, mergedSeq):

	gi1_header=''
	gi2_header=''
	logger.debug("Inside the create merge Seq Obj function: \n\n")
	header1 = seq1Obj.description
	header2 = seq2Obj.description
	logger.debug(header1)
	logger.debug(header2)
	m1 = re.search('gi\|([^\|]+)\|', header1)
	if m1:
		gi1_header = m1.group(1)
	m2 = re.search('gi\|([^\|]+)\|', header2)
	if m2:
		gi2_header = m2.group(1)

	newHeader = header1 + "|ITS-merge|ITS1:" + gi1_header + " ITS2:" + gi2_header
	mergedSeqObj = SeqRecord(Seq(mergedSeq),id = gi1_header, description = newHeader)

	return mergedSeqObj




# append two given sequences (not aligned to each other)
def appendSequences(seq1Obj, seq2Obj):

	appendedSeq = str(seq1Obj.seq + seq2Obj.seq)
	appendedSeqObj = createMergedSeqObj_py(seq1Obj, seq2Obj, appendedSeq)

	return appendedSeqObj


# merge two given aligned sequences (aligned to each other)
def mergeSequences(seq1Obj, seq2Obj):

	indx=0
	mergedSeq_list=[]
	len_seq = len(seq1Obj.seq)
	while indx < len_seq:
		if (seq1Obj.seq[indx] is "-" and seq2Obj.seq[indx] is not "-"):
			mergedSeq_list.append(seq2Obj.seq[indx])
		else:
			mergedSeq_list.append(seq1Obj.seq[indx])
		indx+=1

	merged_seq = str(''.join(mergedSeq_list))
	mergedSeqObj = createMergedSeqObj_py(seq1Obj, seq2Obj, merged_seq)

	return mergedSeqObj

def getChosenSequence(ITS1count, ITS2count, seq1Obj, seq2Obj, append):

	chosenSeq = 'None'

	if (ITS1count > 0 and ITS2count > 0):
		if append:
			chosenSeq = appendSequences(seq1Obj, seq2Obj)
		else:
			chosenSeq = mergeSequences(seq1Obj, seq2Obj)
	elif ITS1count > 0 :
		chosenSeq = seq1Obj
	elif ITS2count > 0 :
		chosenSeq = seq2Obj

	return chosenSeq

def ITS_MSA_py(OutDir, ITS1count, ITS2count, combinedCount,ITS1Fasta, ITS2Fasta, combinedFasta,scriptsDir):

	outputFileName = 'None'
	f_msa_its = open(OutDir+'/ITS_msa.log','w')
	f_msa_its.write("About to adjust direction for its sequences\n")

	# combined sequences exist
	if (combinedCount > 0):
		# Adjusting the direction of the relevant fasta files
		if (ITS1count > 0 or ITS2count > 0):
			os.system("cat %s %s > %s/SEP_ITS1+ITS2.fasta" % (ITS1Fasta,ITS2Fasta,OutDir))
			adjustDirScript = ott_config['general']['OTT_MAIN']+'/ott_scripts/adjustDirFasta.py'
			f_msa_its.write("Calling: python %s -i %s,%s/SEP_ITS1+ITS2.fasta\n" %(adjustDirScript,combinedFasta,OutDir))
			os.system("python %s -i %s,%s/SEP_ITS1+ITS2.fasta" %(adjustDirScript,combinedFasta,OutDir))
		else:
			adjustDirScript = ott_config['general']['OTT_MAIN']+'/ott_scripts/adjustDirFasta.py'
			f_msa_its.write("Calling: python %s -i %s\n" %(adjustDirScript,combinedFasta))
			os.system("python %s -i %s" %(adjustDirScript,combinedFasta))


		# create MSA of the combined sequences (contain both ITS1 and ITS2)
		if (combinedCount > 1):
			combined = OutDir + "/combined.msa"
			#os.system("mafft --auto --ep 0.000000 %s > %s" %(combinedFasta,combined))
			exec_external_command_redirect_output("mafft --auto --ep 0.000000 %s > %s" %(combinedFasta,combined))
		elif(combinedCount == 1):
			combined = OutDir + "/combined.msa"
			shutil.copyfile(combinedFasta,combined)
		else:
			combined = "combined.fasta"

		outputFileName = combined
		# if there are separate sequences, add them to the MSA with 'addfragments' option of MAFFT
		if (ITS1count > 0 or ITS2count > 0):
			os.system("mafft --addfragments %s/SEP_ITS1+ITS2.fasta --multipair %s/combined.msa > %s/combined+sep.msa" %(OutDir,OutDir,OutDir))
			#Check if addfrag worked:
			combined_sep_f = open(OutDir+'/combined+sep.msa','r')
			if is_empty_file(combined_sep_f) != 0:
				outputFileName = OutDir + "/combined+sep.msa"
			else:
				outputFileName = OutDir + "/combined.msa"

	# no combined sequences
	else:
		# MSA for ITS1 sequences (if more than one exist)
		indexIts12FileName  = OutDir + '/ITS_1_2_msa_index.txt'
		with open(indexIts12FileName,'w') as f_index_its:
			f_index_its.write("%s/ITS1_only.msa\n" %OutDir)
			f_index_its.write("%s/ITS2_only.msa\n" %OutDir)
		f_index_its.close()

		if (ITS1count > 1):
			adjustDirScript = ott_config['general']['OTT_MAIN']+'/ott_scripts/adjustDirFasta.py'
			f_msa_its.write("Calling: python %s -i %s\n" %(adjustDirScript,ITS1Fasta))
			os.system("python %s -i %s" %(adjustDirScript,ITS1Fasta))
			os.system("mafft --auto %s > %s/ITS1_only.msa" %(ITS1Fasta,OutDir))
			# system "mafft --retree 2 --maxiterate 0 --bl 62 --op 1.530000 --ep 0.000000 $ITS1Fasta > OutDir/ITS1_only.msa" ; #-bl 62 , --op 1.530000 -> is default, --ep 0.000000 allows large gaps !!!
			outputFileName = OutDir + "/ITS1_only.msa"
		elif(ITS1count == 1):
			outputFileName = ITS1Fasta

		# MSA for ITS2 sequences (if more than one exist)
		if (ITS2count > 1):
			adjustDirScript = ott_config['general']['OTT_MAIN']+'/ott_scripts/adjustDirFasta.py'
			f_msa_its.write("Calling: python %s -i %s\n" %(adjustDirScript,ITS2Fasta))
			os.system("python %s -i %s" %(adjustDirScript,ITS2Fasta))
			#os.system("mafft --auto 0 --ep 0.000000 %s > %s/ITS2_only.msa" %(ITS2Fasta,OutDir))
			os.system("mafft --auto %s > %s/ITS2_only.msa" %(ITS2Fasta,OutDir))

			if outputFileName != 'None':
				os.system("cat %s > %s/SEP_ITS1+ITS2.msa" %(outputFileName,OutDir))
				os.system("cat %s/ITS2_only.msa >> %s/SEP_ITS1+ITS2.msa" %(OutDir,OutDir))

				#Concat ITS1 & ITS2:
				outputFileName = OutDir + "/SEP_ITS1+ITS2.msa"
				prepare_files_for_concat(indexIts12FileName)
				concat_align_command = "perl /groups/pupko/haim/pupkoSVN/trunk/programs/indelReliability/ConcateAlignments.pl %s %s %s NO NA NA" % \
									   (indexIts12FileName, outputFileName, OutDir + "/concat-report.txt")
				os.system(concat_align_command)
				return 'sep_ready'
				#outputFileName = OutDir + "/SEP_ITS1+ITS2.msa"
			else:
				outputFileName = OutDir + "/ITS2_only.msa";
		elif(ITS2count == 1):
			if outputFileName != 'None':
				os.system("cat %s > %s/SEP_ITS1+ITS2.msa" %(outputFileName,OutDir))
				os.system("cat %s >> %s/SEP_ITS1+ITS2.msa" %(ITS2Fasta,OutDir))
				outputFileName = OutDir + "/SEP_ITS1+ITS2.msa"
			else:
				outputFileName = ITS2Fasta


		# Align SEP_ITS1_ITS2: added bu michal to ensure its aligned file
		#os.system("cat %s/SEP_ITS1+ITS2.msa >> %s/SEP_ITS1+ITS2_temp.msa" %(OutDir,OutDir))
		#os.system("mafft --auto --quiet %s/SEP_ITS1+ITS2_temp.msa > %s" %(OutDir,outputFileName))

	return outputFileName

def ITS_CLUSTALO_py(OutDir, ITS1count, ITS2count, combinedCount,ITS1Fasta, ITS2Fasta, combinedFasta,scriptsDir):

	outputFileName = 'None'
	f_msa_its = open(OutDir+'/ITS_msa.log','w')
	f_msa_its.write("About to adjust direction for its sequences\n")

	# combined sequences exist
	if (combinedCount > 0):
		# Adjusting the direction of the relevant fasta files
		if (ITS1count > 0 or ITS2count > 0):
			os.system("cat %s %s > %s/SEP_ITS1+ITS2.fasta" % (ITS1Fasta,ITS2Fasta,OutDir))
			#os.system("cat %s > %s/SEP_ITS1+ITS2.fasta" %(ITS1Fasta,OutDir))
			#os.system("%s >> %s/SEP_ITS1+ITS2.fasta" %(ITS2Fasta,OutDir))

			adjustDirScript = ott_config['general']['OTT_MAIN'] + '/ott_scripts/adjustDirFasta.py'
			f_msa_its.write("Calling: python %s -i %s,%s/SEP_ITS1+ITS2.fasta\n" % (adjustDirScript, combinedFasta, OutDir))
			os.system("python %s -i %s,%s/SEP_ITS1+ITS2.fasta" % (adjustDirScript, combinedFasta, OutDir))
		else:
			adjustDirScript = ott_config['general']['OTT_MAIN'] + '/ott_scripts/adjustDirFasta.py'
			f_msa_its.write("Calling: python %s -i %s\n" % (adjustDirScript, combinedFasta))
			os.system("python %s -i %s" % (adjustDirScript, combinedFasta))

		# create MSA of the combined sequences (contain both ITS1 and ITS2) - ClastaO
		if (combinedCount > 1):
			combined = OutDir + "/combined.msa"
			MSA_cmd = '%s -i %s -o %s --outfmt=fasta' % (ott_config['diff_soft']['ClastaO'], combinedFasta, combined)
			exec_external_command_redirect_output(MSA_cmd)
		elif(combinedCount == 1):
			combined = OutDir + "/combined.msa"
			shutil.copyfile(combinedFasta,combined)
		else:
			combined = "combined.fasta"

		outputFileName = combined
		# if there are separate sequences, add them to the MSA with 'addfragments' option of MAFFT
		if (ITS1count > 0 or ITS2count > 0):
			os.system("mafft --addfragments %s/SEP_ITS1+ITS2.fasta --multipair %s/combined.msa > %s/combined+sep.msa" %(OutDir,OutDir,OutDir))
			#Check if addfrag worked:
			combined_sep_f = open(OutDir+'/combined+sep.msa','r')
			if is_empty_file(combined_sep_f) != 0:
				outputFileName = OutDir + "/combined+sep.msa"
			else:
				outputFileName = OutDir + "/combined.msa"

	# no combined sequences
	else:
		# MSA for ITS1 sequences (if more than one exist)
		if (ITS1count > 1):
			adjustDirScript = ott_config['general']['OTT_MAIN']+'/ott_scripts/adjustDirFasta.py'
			f_msa_its.write("Calling: python %s -i %s\n" %(adjustDirScript,ITS1Fasta))
			os.system("python %s -i %s" %(adjustDirScript,ITS1Fasta))
			os.system("mafft --auto %s > %s/ITS1_only.msa" %(ITS1Fasta,OutDir))
			# system "mafft --retree 2 --maxiterate 0 --bl 62 --op 1.530000 --ep 0.000000 $ITS1Fasta > OutDir/ITS1_only.msa" ; #-bl 62 , --op 1.530000 -> is default, --ep 0.000000 allows large gaps !!!
			outputFileName = OutDir + "/ITS1_only.msa"
		elif(ITS1count == 1):
			outputFileName = ITS1Fasta

		# MSA for ITS2 sequences (if more than one exist)
		if (ITS2count > 1):
			adjustDirScript = ott_config['general']['OTT_MAIN']+'/ott_scripts/adjustDirFasta.py'
			f_msa_its.write("Calling: python %s -i %s\n" %(adjustDirScript,ITS2Fasta))
			os.system("python %s -i %s" %(adjustDirScript,ITS2Fasta))
			os.system("mafft --auto %s > %s/ITS2_only.msa" %(ITS2Fasta,OutDir))

			if outputFileName:
				os.system("cat %s > %s/SEP_ITS1+ITS2.msa" %(outputFileName,OutDir))
				os.system("cat %s/ITS2_only.msa >> %s/SEP_ITS1+ITS2.msa" %(OutDir,OutDir))
				outputFileName = OutDir + "/SEP_ITS1+ITS2.msa"
			else:
				outputFileName = OutDir + "/ITS2_only.msa";
		elif(ITS2count == 1):
			if outputFileName:
				os.system("cat %s > %s/SEP_ITS1+ITS2.msa" %(outputFileName,OutDir))
				os.system("cat %s >> %s/SEP_ITS1+ITS2.msa" %(ITS2Fasta,OutDir))
				outputFileName = OutDir + "/SEP_ITS1+ITS2.msa"
			else:
				outputFileName = ITS2Fasta


		# Align SEP_ITS1_ITS2: added bu michal to ensure its aligned file
		os.system("cat %s/SEP_ITS1+ITS2.msa >> %s/SEP_ITS1+ITS2_temp.msa" %(OutDir,OutDir))
		os.system("mafft --auto --quiet %s/SEP_ITS1+ITS2_temp.msa > %s" %(OutDir,outputFileName))

	return outputFileName



#Perform Blast on fasta (with one ITS type) to choose the representetive seq of this type.
# The other seqs will be removed:
def pickSeqRecord_py(context,fileName,outDir):

	#Need to add the code  so we'll have the length of each type in case the blast results are not definit we
	# can use this criteria to choose the representetiv seq
	accession_id_list=[] # list of accessions to be blasted
	records_list = list(SeqIO.parse(fileName, "fasta"))
	missing_accessions = dict()  # 0 missing, 1 exist
	accession_id_vs_SeqRecord = dict()
	for seq in records_list:
		acc_id = getPropertyFromFastaSeqHeader(seq.description, "gi")
		accession_id_list.append(acc_id)
		accession_id_vs_SeqRecord[acc_id] = seq
	logger.debug("accession_id_list")
	logger.debug(accession_id_list)

	if len(records_list) > 1:
		logger.debug("File name - %s, len(records_list) - %d" %(fileName,len(records_list)))
		type_name = os.path.basename(fileName) # combined/ITS1/ITS2
		blast_file = outDir + '/' + type_name + '_all_v_all.blastn'
		os.system("formatdb -i %s -pF" %fileName)
		os.system("blastall -p blastn -d %s -i %s -v 100000 -b 100000 -e 1e-5 -m 8 > %s" %(fileName,fileName,blast_file))
		f_blastResults = open(blast_file)
		logger.debug("Generating dictionary of bitscores between seqs according to %s" % blast_file)

		dr = csv.DictReader(f_blastResults, delimiter='\t',fieldnames=['query_id','subject_id','pct_identity','align_len','mismatches',
															  'gap_openings','q_start','q_end','s_start','s_end','eval','bitscore'])
		rowslist = list(dr)

		#For each pair save the score in dict: check that all accession are here:
		accession_pair_vs_MaxScore=dict()
		pairs_list=[]
		min_score=999999
		for row in rowslist:
			score_list=[]
			query_id=get_acc_from_id(row['query_id'])
			subj_id = get_acc_from_id(row['subject_id'])
			pair_key = get_taxon_gis_key(query_id,subj_id)
			pairs_list.append(pair_key)
			score = float(row['bitscore'])
			if score < min_score:	#save minimum score:
				min_score = score
			#save all score results of this pair as a list:
			if pair_key not in accession_pair_vs_MaxScore:
				score_list.append(float(score))
				accession_pair_vs_MaxScore[pair_key] = score_list
			else:
				score_list = list(accession_pair_vs_MaxScore[pair_key])
				score_list.append(float(score))
				accession_pair_vs_MaxScore[pair_key] = score_list
			#logger.debug(score_list)
		#logger.debug("accession_pair_vs_MaxScore:")
		#for key in accession_pair_vs_MaxScore:
		#	logger.debug(key)
		#	logger.debug(accession_pair_vs_MaxScore[key])
		logger.debug("pairs_list:")
		logger.debug(pairs_list)
		logger.debug("min score:")
		logger.debug(min_score)

		#calculating the max, of all avg values. The value is an avg of all scores of each accession.
		#In case there are no results for some pairs we need to set them as min score
		max_avg_accessionId=0
		max_avg=0
		for accession_id in accession_id_list:
			accession_list_withoutHead = list(accession_id_list)
			accession_list_withoutHead.remove(accession_id)
			sub_id_list=[]
			logger.debug("accession_list_withoutHead:")
			logger.debug(accession_list_withoutHead)
			logger.debug("Calc avg for accession id: %s" %accession_id)
			for pair in pairs_list:
				query_id = pair.split('$')[0]
				sub_id = pair.split('$')[1]
				if accession_id in query_id and accession_id not in sub_id:
					sub_id_list.insert(0,sub_id)
					score_list = list(accession_pair_vs_MaxScore[pair])
					avg_total=0;cnt=0
					for score in score_list:
						avg_total+=score
						cnt+=1
					if avg_total != 0:
						avg=avg_total/cnt
					if avg > max_avg:
						max_avg = avg
						max_avg_accessionId=accession_id
				else: # Only result for the accession vs itself
					score_list = list(accession_pair_vs_MaxScore[pair])
					avg_total=0;cnt=0
					for score in score_list:
						avg_total+=score
						cnt+=1
					if avg_total != 0:
						avg=avg_total/cnt
					if avg > max_avg:
						max_avg = avg#
						max_avg_accessionId=accession_id

		logger.debug("max_avg, accession_id")
		logger.debug("%s,%s" %(max_avg,max_avg_accessionId))

		logger.debug("Selected Accession with highest blast Bitscore is %s" % max_avg_accessionId)
		#check if max_avg_accessionId is 0 -> then we need to select according to length:
		#This case happens when there are no blast results for any of the accessions pairs:
		if max_avg_accessionId == 0:
			accession_id_longest = return_longest_seq(context,accession_id_list)
			logger.debug("(length criteria) Found selected sequence with accession id %s" % accession_id_longest)
			return accession_id_vs_SeqRecord[accession_id_longest]
		else:
			for seq_record in SeqIO.parse(fileName, "fasta"):
				accession_id = getPropertyFromFastaSeqHeader(seq_record.description, "gi")
				logger.debug(accession_id)
				logger.debug(max_avg_accessionId)
				if str(accession_id) == str(max_avg_accessionId):
					selectedSeq = seq_record
					logger.debug("Found selected sequence with accession id %s" % accession_id)
					return selectedSeq
	else:
		return 0

def getTheFirstSequence_py(fileName):

	for seq_record in SeqIO.parse(fileName, "fasta"):
		return seq_record

def pickFromFile_py(context,count,fileName,outDir):

	if (count > 1):
		seq_record = pickSeqRecord_py(context,fileName,outDir)
	elif(count == 1):
		seq_record = getTheFirstSequence_py(fileName)

	logger.debug(seq_record)

	if seq_record is 0:
		logger.debug("seq_record is 0")
	else:
		logger.debug("seq_record %s" %seq_record.description)
	return seq_record


def pickOneSeqPerSpeciesMSA(context, inputFile, outDir, outputFile, scriptsDir, logFile,append):

	f_log = open(logFile,'a')
	f_log.write ("Calling2 perl %s/SplitFastaSeqsBySpecies.pl %s %s/species\n" %(scriptsDir,inputFile,outDir));
	os.system("perl %s/SplitFastaSeqsBySpecies.pl %s %s/species" %(scriptsDir,inputFile,outDir))

	f_log.write("context.its_accession_vs_type Dictionary:")
	for key in context.its_accession_vs_type.keys():
		f_log.write('%s: %s\n' %(key,context.its_accession_vs_type[key]))

	f_log.write("outputFile:\n")
	f_log.write(outputFile)
	f_log.write('\n')
	Final_records=[]
	for taxon_fasta_file in glob.glob(outDir+"/species/*.fasta"):

		ITS1count = 0
		ITS2count = 0
		combinedCount = 0

		base_file_name = os.path.basename(taxon_fasta_file)
		speciesDirName = base_file_name.replace('.fasta', '')
		speciesDir = outDir+"/species/" +speciesDirName

		seq_records = SeqIO.parse(taxon_fasta_file,'fasta')
		for seq in seq_records:
			f_log.write(seq.id)
			seq_accession = getPropertyFromFastaSeqHeader(seq.description, "gi")
			if seq_accession in context.its_accession_vs_type:
				f_log.write(context.its_accession_vs_type[seq_accession])
				if context.its_accession_vs_type[seq_accession] is 'combined':
					SeqIO.write(seq, speciesDir + "/combined", "fasta")
					combinedCount+=1
				if context.its_accession_vs_type[seq_accession] is 'its1':
					SeqIO.write(seq, speciesDir + "/ITS1_only", "fasta")
					ITS1count+=1
				if context.its_accession_vs_type[seq_accession] is 'its2':
					SeqIO.write(seq, speciesDir + "/ITS2_only", "fasta")
					ITS2count+=1

		#Check the flow was correct and there is only one seq per type
		if (ITS1count > 1 or ITS2count > 1 or combinedCount > 1):
			f_log.write("More than a single ITS per type was found in %s" %taxon_fasta_file)
			return Fail
		if ((ITS1count > 1 or ITS2count > 1) and (combinedCount > 0)):
			f_log.write("Combined ITS was found as well as ITS1/ITS2 in %s" % taxon_fasta_file)
			return Fail

		f_log.write("taxon_fasta_file = %s)\n" % taxon_fasta_file)
		f_log.write("Counters are: %d,%d,%d\n" %(combinedCount,ITS1count,ITS2count))
		if combinedCount > 0 :
			# pick a combined sequence if available:
			chosenSeq = pickFromFile_py(context,combinedCount, speciesDir +"/combined", speciesDir)
			Final_records.append(chosenSeq)
			#SeqIO.write(chosenSeq, outputFile, "fasta")
		elif (ITS1count > 0 and ITS2count > 0):
			# else, merge two separated sequences - ITS1 + ITS2:
			seqObj1 = getTheFirstSequence_py(speciesDir +"/ITS1_only")
			seqObj2 = getTheFirstSequence_py(speciesDir +"/ITS2_only")
			f_log.write("Merging/appending ITS1 and ITS2\n")
			chosenSeq = getChosenSequence(1, 1, seqObj1, seqObj2, append)
			Final_records.append(chosenSeq)
			f_log.write("species: ITS1 and ITS2 were merged\n")
		elif (ITS1count > 0) :
			chosenSeq = getTheFirstSequence_py(speciesDir + "/ITS1_only")
			Final_records.append(chosenSeq)
			f_log.write("Using ITS1\n")
		elif (ITS2count > 0) :
			chosenSeq = getTheFirstSequence_py(speciesDir + "/ITS2_only")
			Final_records.append(chosenSeq)
			f_log.write("Using ITS2\n")
		else:
			f_log.write("ERROR - species has no ITS sequence\n")

	#Write all chosen sequences to the final msa:
	SeqIO.write(Final_records, outputFile, "fasta")

	return



#This function split ITS sequences per species, per type and choose a representative sequence for each type:
def pickOneITSTypePerSpeciesFasta_py(context,inputFile, outDir, outputFile, scriptsDir, logFile):

	f_log = open (logFile,'w')
	#Also save the seqs to seperate fasta files:
	total_combined_f = open(outDir+'/combined.fasta','w')
	total_ITS1_f = open(outDir+'/ITS1_only.fasta','w')
	total_ITS2_f = open(outDir+'/ITS2_only.fasta','w')
	total_combined_cnt = 0
	total_its1_cnt = 0
	total_its2_cnt = 0

	f_log.write("Calling1 perl %s/SplitFastaSeqsBySpecies.pl %s %s/species_all\n" %(scriptsDir,inputFile,outDir))
	os.system("perl %s/SplitFastaSeqsBySpecies.pl %s %s/species_all" %(scriptsDir,inputFile,outDir))

	chosenSequences = []

	speciesFiles = glob.glob(outDir+"/species_all/*.*")    #   grep { $_ ne '.' && $_ ne '..' } readdir(SP);
	#out_rec = seq_sep_records = SeqIO.parse(outputFile, 'fasta')    #Bio::SeqIO->new("-file" => ">$outputFile", "-format" => "Fasta");
	speciesFiles_temp = [w.replace(outDir+"/species_all/", '') for w in speciesFiles]
	TaxIDs_ITS_list = [w.replace(".fasta", '') for w in speciesFiles_temp]
	logger.debug(TaxIDs_ITS_list)

	species_counter = 0
	chosenSequences = []
	f_log.write("Processing $species\n")
	for seq_file in speciesFiles:
		if('.fasta' in seq_file):
			base_file_name = os.path.basename(seq_file)
			speciesDirName = base_file_name.replace('.fasta','')
		#Find Accessions:
		#records = SeqIO.parse(seq_file,'fasta')
		#for seq in records:
		#	accession_id = getPropertyFromFastaSeqHeader(seq.description, "gi")
		#	context.its_accession_ids.append(accession_id)
		#logger.debug("context.its_accession_ids")
		#logger.debug(context.its_accession_ids)
		speciesDir = "%s/species_all/%s" %(outDir,speciesDirName)
		speciesForMsaDir = "%s/species/%s" %(outDir,speciesDirName) # For later use

		create_dir_if_not_exists(speciesDir)
		create_dir_if_not_exists(speciesForMsaDir)

		(ITS1count, ITS2count, combinedCount) = splitITS_py(speciesDirName,context, speciesDir+'.fasta', speciesDir + "/ITS1_only",speciesDir + "/ITS2_only", speciesDir + "/combined")
		f_log.write("Completed 1st splitITS_py: ITS1count=%d, ITS2count=%d, combinedCount =%d)\n" % (ITS1count, ITS2count, combinedCount))
		context.its_taxa_vs_counts[speciesDir] = [ITS1count, ITS2count, combinedCount]

	#calc longest in case of no blast results:
	calc_longest_accordingToType(context)

	#For each species select representative seq:
	for speciesDir in context.its_taxa_vs_counts:
		ITS1count = context.its_taxa_vs_counts[speciesDir][0]
		ITS2count = context.its_taxa_vs_counts[speciesDir][1]
		combinedCount = context.its_taxa_vs_counts[speciesDir][2]
		if (combinedCount > 0):
			chosen_seq = pickFromFile_py(context,combinedCount, speciesDir + "/combined", speciesDir)
			if chosen_seq is not 0:
				species_counter+=1
				chosenSequences.append(chosen_seq)
				SeqIO.write(chosen_seq, total_combined_f, "fasta")
				total_combined_cnt+=1

		# Else,  merge two separated sequences - ITS1 + ITS2:
		elif(ITS1count > 0 or ITS2count > 0):

			if (ITS1count > 0):
				chosen_seq = pickFromFile_py(context,ITS1count, speciesDir + "/ITS1_only", speciesDir)
				if chosen_seq is not 0:
					species_counter += 1
					chosenSequences.append(chosen_seq)
					SeqIO.write(chosen_seq, total_ITS1_f, "fasta")
					total_its1_cnt += 1

			if (ITS2count > 0):
				chosen_seq = pickFromFile_py(context,ITS2count, speciesDir + "/ITS2_only", speciesDir)
				if chosen_seq is not 0:
					species_counter += 1
					chosenSequences.append(chosen_seq)
					SeqIO.write(chosen_seq, total_ITS2_f, "fasta")
					total_its2_cnt += 1


	f_log.write("A total of %d sequences were chosen\n" %len(chosenSequences))
	f_log.write("Number of species =%d\n" %species_counter)
	f_log.write("Total seq count according to ITS type: combined=%d, its1=%d, its2=%d\n" %(total_combined_cnt,total_its1_cnt,total_its2_cnt))
	f_out = open(outputFile, 'w')
	for seq in chosenSequences:
		SeqIO.write(seq, f_out, "fasta")

	total_combined_f.close()
	total_ITS1_f.close()
	total_ITS2_f.close()

	return (total_its1_cnt,total_its2_cnt,total_combined_cnt)


def formatKeyWords(features_concat):

	logger.debug(features_concat)

	internal_transcribed = re.compile("internal transcribed spacer ", re.IGNORECASE)
	internal_trasncribed = re.compile("internal trasncribed spacer ", re.IGNORECASE)
	its_1 = re.compile("its 1", re.IGNORECASE)
	its_2 = re.compile("its 2", re.IGNORECASE)
	its1 = re.compile(" its1", re.IGNORECASE)
	its2 = re.compile(" its2", re.IGNORECASE)
	its = re.compile(" its ", re.IGNORECASE)

	features_concat_new = internal_transcribed.sub("ITS", features_concat)
	features_concat_new = internal_trasncribed.sub("ITS", features_concat_new)
	features_concat_new = its_1.sub("ITS1", features_concat_new)
	features_concat_new = its_2.sub("ITS2", features_concat_new)
	features_concat_new = its1.sub(" ITS1", features_concat_new)
	features_concat_new = its2.sub(" ITS2", features_concat_new)
	features_concat_new = its.sub(" ITS ", features_concat_new)

	logger.debug(features_concat_new)

	return features_concat_new


def splitITS_py(taxId, context,input, its1Output, its2Output, combinedOutput):

	its1Out = open (its1Output,'w')
	its2Out = open (its2Output,'w')
	combinedOut = open (combinedOutput,'w')
	logger.debug("Path for combined file:")
	logger.debug(combinedOutput)

	ITS1count=0
	ITS2count=0
	combinedCount=0

	seq_count=0
	#Find accession numbers of ITS seqs

	#seq_count = context.taxon_Id_df.index(taxId)
	idx_list=[]
	all_idx_taxId = context.taxon_Id_df[context.taxon_Id_df == taxId]
	accession_indxes = all_idx_taxId.index.tolist()

	#while seq_count < len(context.Accession_df):
	for seq_count in accession_indxes:
		if str(context.Accession_df[seq_count]) in context.its_accession_ids:
			taxon_fasta_str = '/' + str(context.taxon_Id_df[seq_count]) + '.fasta' # Added / to solve BUG_md#1
			if taxon_fasta_str in input:
				features_concat = context.definition_df[seq_count] + '|' + context.gene_name_location_df[seq_count] + '|' + \
								  context.note_df[seq_count] + '|' + context.mol_type_df[seq_count] + '|' + context.note_df[seq_count] + \
								  '|' + context.product_df[seq_count]
				desc = formatKeyWords(features_concat)

				#Need to get all the features data from features feild
				temp_its_accession = str(context.Accession_df[seq_count])
				if ('ITS1' in desc) and ('ITS2' in desc):
					combinedOut.write(context.desc_df[seq_count]+'\n')
					combinedOut.write(context.sequenceData_df[seq_count]+'\n')
					seq_len = len(context.sequenceData_df[seq_count])
					combinedCount+=1
					#Save Accession ITS type
					context.its_accession_vs_type[temp_its_accession] = 'combined'
					context.its_accession_vs_length[temp_its_accession] = seq_len

				elif ('ITS1' in desc) and ('ITS2' not in desc):
					its1Out.write(context.desc_df[seq_count] + '\n')
					its1Out.write(context.sequenceData_df[seq_count] + '\n')
					seq_len = len(context.sequenceData_df[seq_count])
					ITS1count += 1
					#Save Accession ITS type
					context.its_accession_vs_type[temp_its_accession] = 'its1'
					context.its_accession_vs_length[temp_its_accession] = seq_len

				elif ('ITS1' not in desc) and ('ITS2' in desc):
					its2Out.write(context.desc_df[seq_count] + '\n')
					its2Out.write(context.sequenceData_df[seq_count] + '\n')
					seq_len = len(context.sequenceData_df[seq_count])
					ITS2count += 1
					# Save Accession ITS type
					context.its_accession_vs_type[temp_its_accession] = 'its2'
					context.its_accession_vs_length[temp_its_accession]=seq_len

	return (ITS1count, ITS2count, combinedCount)


#-------------------------------------------------------------------------------------------------------
#
#							MAIN  - ITS flow (converted from perl)
#
#-------------------------------------------------------------------------------------------------------
def main_ITS_py(context,OutDir, fasta, gene_db, scriptsDir, configFile, guidanceFlag, msa_software):

	f_mainITS_log = open(OutDir+'/MainITS.log','w')
	pickFromFastalog = OutDir+"/pickOneSeqFromFasta.log"
	pickFromMSAlog = OutDir+"/pickOneSeqFromMSA.log"
	its1fasta = OutDir+"/ITS1_only.fasta"
	its2fasta = OutDir+"/ITS2_only.fasta"
	itscombfasta = OutDir+"/combined.fasta"

	create_dir_if_not_exists(OutDir)

	#Create list of Accessions of ITS:
	logger.debug("Create list of AccessionId for all ITS sequences - start")
	records_all = SeqIO.parse(fasta, 'fasta')
	for seq in records_all:
		accession_id = getPropertyFromFastaSeqHeader(seq.description, "gi")
		context.its_accession_ids.append(accession_id)
	logger.debug("Number of AccessionId for all ITS sequences - %d" %len(context.its_accession_ids))
	logger.debug(context.its_accession_ids)

	f_mainITS_log.write("Filtering %s - each species will have at most one ITS seq of each type (ITS1/ITS2/combined) \n" %fasta)
	#--> Need to convert to python -> removed gbIndex since we use Database instead
	(ITS1count, ITS2count, combinedCount) = pickOneITSTypePerSpeciesFasta_py(context,fasta, OutDir, OutDir+"/oneITSTypePerSpecies.fasta", scriptsDir, pickFromFastalog)
	f_mainITS_log.write("ITS1count=%d, ITS2count=%d, combinedCount =%d, )" % (ITS1count, ITS2count, combinedCount))
	f_mainITS_log.write("Completed - pickOneITSTypePerSpeciesFasta_py\n")

	append = 1
	if combinedCount > 0 :
		append = 0   # merge
	f_mainITS_log.write("append flag = %d (if 0-> merge,if 1-> append)" %append)

	if(ITS1count > 0 or ITS2count > 0 or combinedCount > 0):

		# Do MSA for the ITS seqs
		f_mainITS_log.write("Found ITS seqs - starting MSA for ITS\n")


		if (msa_software == 'ClustalOmega'):
			f_mainITS_log.write("MSA software -> ClustalOmega\n")
			MSAfileName = ITS_CLUSTALO_py(OutDir, ITS1count, ITS2count, combinedCount,its1fasta, its2fasta, itscombfasta,scriptsDir)
		else:
			f_mainITS_log.write("MSA software -> MAFFT\n")
			MSAfileName = ITS_MSA_py(OutDir, ITS1count, ITS2count, combinedCount,its1fasta, its2fasta, itscombfasta,scriptsDir)

		if MSAfileName is 'sep_ready':
			logger.debug("ITS from ITS1 and ITS2, no combined data")
			shutil.copyfile(OutDir + '/SEP_ITS1+ITS2.msa', OutDir+"/oneSeqPerSpecies.msa")
			return
		f_mainITS_log.write("After MSA - selecting a single ITS seq per species. MSA file: %s\n" % MSAfileName)
		pickOneSeqPerSpeciesMSA(context,MSAfileName, OutDir, OutDir+"/oneSeqPerSpecies.msa", scriptsDir, pickFromMSAlog,append)

		logger.debug(os.path.exists(OutDir + "/combined+sep.msa"))
		if os.path.exists(OutDir + "/combined+sep.msa"):
			combined_sep_f = open(OutDir + '/combined+sep.msa', 'r')
			if is_empty_file(combined_sep_f) != 0:
				sep_fasta_f = ("%s/SEP_ITS1+ITS2.fasta" %OutDir)
				logger.debug("About to count number of seqs in %s\n" %sep_fasta_f)

				logger.debug("GUIDANCE Flag is set to: %s\n" %guidanceFlag)

				if (guidanceFlag == "GUIDANCE"):
					logger.debug("GUIDANCE is running on ITS\n")
					guidance_for_its(context,OutDir + "/combined.fasta", sep_fasta_f, OutDir)
				elif (guidanceFlag == "Trimal"):
					Trimal_cf = context.UserFlags_dict['Trimal_CutOff']
					runTrimal(OutDir + "/combined+sep.msa", Trimal_cf)
					logger.info("Files after Trimal at: %s" % OutDir + "/combined+sep.msa")
				else:
					logger.debug("Guidance Flag is set to False !!!")
			elif os.path.exists(OutDir + "/combined.msa"):
				if (guidanceFlag == "GUIDANCE"):
					logger.debug("GUIDANCE is running on ITS\n")
					guidance_for_its(context, OutDir + "/combined.fasta", 'None', OutDir)
				elif (guidanceFlag == "Trimal"):
					Trimal_cf = context.UserFlags_dict['Trimal_CutOff']
					runTrimal(OutDir + "/combined.msa", Trimal_cf)
					logger.info("Files after Trimal at: %s" % OutDir + "/combined.msa")
				else:
					#TO DO - execute GUIDACE even in this case - it should be standatd GUIDANCE execution and not special like the above
					logger.debug("NOT Calling GUIDANCE since combined file doesn't exists\n")

	return

