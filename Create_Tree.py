from ploidbCommon import *
import time
import csv
import sys
import ntpath
import os
import re
import mmap
import argparse
from Bio import AlignIO
import fileinput


import run_jModelTest
import createJobFile
from parse_jModelTest2_file import *
from MBconfig_Dictionary import mb_block_dict

PRE_TREE_DIR = 'PreTree_Dir'
PHY_EXT = '_phy.phy'
JMT_EXT = '_phy.phy.jmt'
FASTA_FORMAT = "fasta"
PHYLIP_FORMAT = "phylip-relaxed"

TABLE_HEADERS = "Model\s+-lnL\s+K\s+TEST\s+delta\s+weight\s+cumWeight(\s+uDelta)?"
TREE_LINE = "Tree for the best TEST model = "

DEBUG_PRINT = 1


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Run external process
#
def exec_external_command_redirect_output_Tree(command_to_exec, log_file, outfile=None, errfile=None, cwd=None):
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

	log_file.write("Executing the following command %s - out written to %s error written to %s", command_to_exec, outfile, errfile)
	log_file.write("cwd=%s" % cwd)
	p = subprocess.Popen(command_to_exec, shell=True, cwd=cwd, stdout=command_out, stderr=command_err)
	# Added by Michal - Stuck in some cases (large input files for formatdb
	stdout, stderr = p.communicate()
	retval = p.wait()

	if retval == 0:
		log_file.write("Execution return code was %i for command %s" % (retval, command_to_exec))
	else:
		log_file.write("Execution return code was %i for command %s" % (retval, command_to_exec))

	if out_handle is not None:
		out_handle.close()
	if err_handle is not None:
		err_handle.close()

	return retval

# CDHITcommand_to_exec = 'cd-hit -i %s -d 0 -o %s -c 0.98 -n 5 -G 1 -g 1 -b 25 -s 0.0 -aL 0.0 -aS 0.0 -T 4 -M 16000' % (context.fasta_cdhit_temp_input_f,context.fasta_cdhit_temp_output_f)
# exec_external_command_redirect_output(command_to_exec=CDHITcommand_to_exec,outfile=context.working_dir + "/CD-HIT_perTaxID.out",errfile=context.working_dir + "/CD-HIT_perTaxID.err")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def convert_fasta_to_phylip(input_file, output_file):
	with open(input_file, "rU") as input_handle:
		alignments = AlignIO.parse(input_handle, FASTA_FORMAT)
		with open(output_file, "w") as output_handle:
			AlignIO.write(alignments, output_handle, PHYLIP_FORMAT)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def get_Nparams_ForModel(Info_Creteria, file_reader):
	# Info_Creteria = information criteria : AIC AICc...
	#Models: GTR+I etc...
	print(Info_Creteria)
	print(file_reader)

	#extract tree features
	test_tree_line = re.sub("TEST", Info_Creteria, TREE_LINE)
	if re.search("(?="+test_tree_line+").+", file_reader) is None:
		return 11
	else:
		re.search("(?="+test_tree_line+").+", file_reader).group(0)
		selected_model_tree = re.search("(?="+test_tree_line+").+", file_reader).group(0)
	#test_tree_features = parse_newick_tree.parse_tree_string(selected_model_tree)              # mdBUG Check the meaning of this line - Is needed??????

	#extract features from table
	test_table_headers = re.sub("TEST", Info_Creteria, TABLE_HEADERS)
	test_headers_match = re.search(test_table_headers + "\r?\n\-+\s\r?\n", file_reader, re.M)
	if test_headers_match is None:
		return None
	test_table_1st_line = re.search(".*", file_reader[test_headers_match.end():]).group(0)
	selected_model_line_features = re.split("\s+", test_table_1st_line)
	# The K column
	Info_Creteria_Nparams = selected_model_line_features[2]
	Sel_Model=selected_model_line_features[0]
	#selected_model = re.sub("\+", "", selected_model) #remove + from model name (.+I+G ==> .IG)
	return (Info_Creteria_Nparams,Sel_Model)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def fixed_model_and_create_blk(phy_file,cluster_dir,User_Model):

	#User_Model = "GTR+G"
	logger.debug("--------------------------------------------------")
	logger.debug("         Selected Model is (User Model): %s" %User_Model)
	logger.debug("--------------------------------------------------")

	infile = open(phy_file, 'r')
	firstLine = infile.readline()
	data=firstLine.split(' ')
	data[-1] = data[-1].strip() # remove \n
	# Length of Alignment:
	Nchars=int(data[2])          # Alignment length
	logger.debug("Nchars (length of Alignment) is %d" %Nchars)

	# Create the MB data blockes from the jModel output files:
	py_MBblk_f = open(phy_file + ".jmt_mb_BLK.txt",'w')
	py_MBblk_f.write("\nBEGIN MRBAYES;\n")
	py_MBblk_f.write(" Charset %s = 1 - %s;\n" % (phy_file, str(Nchars)))
	py_MBblk_f.write(" Partition Dummy = 1:%s;\n" % phy_file)
	py_MBblk_f.write(" Set partition = Dummy;\n")
	py_MBblk_f.write(mb_block_dict(User_Model) + "\n")
	py_MBblk_f.write("END;\n\n\n")

	return 0

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def check_model_and_create_blk(phy_file,jmt_file,cluster_dir):

	#Get the length of the alignment (Nchars) from the Phy file
	status = 0
	infile = open(phy_file, 'r')
	firstLine = infile.readline()
	data=firstLine.split(' ')
	data[-1] = data[-1].strip() # remove \n
	# Length of Alignment:
	Nchars=int(data[2])          # Alignment length
	print("Nchars !!!!!!!!!!!!!! %d" %Nchars)
	# Get the number of parameters per model (Nparams) from the jModelTest output file:
	# First calc (Nchar/Nparams_AIC < 40) for AIC, if false choose AICc
	#----------------------------------------------------------------------------------
	jmtfile = open(jmt_file, 'r')
	file_reader = jmtfile.read()
	print("MMMMMMMMMMMMMMMMMMMMMMMM")
	#print('AIC',file_reader[0])
	if get_Nparams_ForModel('AIC',file_reader) == 11:
		return 1
	Nparams_AIC=int(get_Nparams_ForModel('AIC',file_reader)[0])
	print("Nchars=", Nchars, ", Nparams=", Nparams_AIC)
	if (Nchars/Nparams_AIC <= 40):
		Model='AIC'
	else:
		Model='AICc'
	print("--------------------------------------------------")
	print("         Selected Model is: ", Model)
	print("--------------------------------------------------")
	# Add the parsing for looking for the right model in Model:
	if get_Nparams_ForModel(Model,file_reader) == 11:
		return 1
	Selected_Model = get_Nparams_ForModel(Model,file_reader)[1]
	#print("###########$$$$$$$$$$$$$$$$$$$$ %s" % Selected_Model)
	#mrB_outDir = cluster_dir
	# Create the MB data blockes from the jModel output files:
	#command_MrBayesBlk = "perl " + SCRIPTS_DIR + "MrBayesBlk_fromjMT.pl " + jmt_file + " " + Selected_Model + " " + str(Nchars)
	#os.system(command_MrBayesBlk)

	# Create the MB data blockes from the jModel output files:
	py_MBblk_f = open(jmt_file + "_mb_BLK.txt",'w')
	py_MBblk_f.write("\nBEGIN MRBAYES;\n")
	py_MBblk_f.write(" Charset %s = 1 - %s;\n" % (jmt_file, str(Nchars)))
	py_MBblk_f.write(" Partition Dummy = 1:%s;\n" % jmt_file)
	py_MBblk_f.write(" Set partition = Dummy;\n")
	py_MBblk_f.write(mb_block_dict(Selected_Model) + "\n")
	py_MBblk_f.write("END;\n\n\n")

	return 0

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Check which model to use according to jMT file:
def part2_create_mrbayes_blk(genus, fastaFilesList):

	with open(fastaFilesList, "r") as pathFile:
		for FastaPath in pathFile:
			file_no_extension, file_extension = os.path.splitext(FastaPath)	# remove file extension
			#fasta_file=FastaPath.replace('\n','')							# Remove \n from path

			phy_file=file_no_extension+PHY_EXT
			jmt_file_path=phy_file + ".jmt"			# jmt output file path + file name
			if(DEBUG_PRINT==1): print(phy_file)
			if(DEBUG_PRINT==1): print(jmt_file_path)

			# Cluster dir path - folder only
			cluster_dir = os.path.dirname(FastaPath)
			if(DEBUG_PRINT==1): print(cluster_dir)
			if (check_model_and_create_blk(phy_file,jmt_file_path,(cluster_dir+'/'))) == 1:
				return 1

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def create_mb_config_file(genus_name,working_dir,fastaFilesList,Ngen_param,scripts_dir,clock_model,relaxed_branch,split_flag):
	#################################################################
	# Create cluster_to_concat_list.txt file and mb config file
	#################################################################
	try:
		concatDir = (working_dir + "/concat/")
		os.chdir(concatDir)
		f_w = open("cluster_to_concat_list.txt", 'w')
		f_r = open("fasta-files-to-concat.txt", 'r')
	except ValueError:
		print("Can't create file cluster_to_concat_list.txt at %s" %concatDir)
	if not os.path.isfile("fasta-files-to-concat.txt"):
		concat_aliged_file=cluster_dir + "/seqs-organism-concat.fasta "
	else:
		with open("fasta-files-to-concat.txt", 'r') as f_r:
			outgroup_seq_file=0;
			num_of_concat_clusters=0
			for line in f_r:
				NoNewLine = line.replace("\n", "")
				f_w.write(NoNewLine + ";")
				# Check for outgroup files:
				# Check for chloroplast
				if "chloroplast-all" in NoNewLine:
					NoNewLine = line.replace("concat/", "999/")
					cluster_dir = os.path.dirname(NoNewLine)
					if(os.path.isfile(cluster_dir + "/outseq-organism.fasta")):
						f_w.write(cluster_dir + "/outseq-organism.fasta;")
						if outgroup_seq_file==0:
							Nex_outgroupPath=(cluster_dir + "/outseq-organism.fasta")
							outgroup_seq_file=1
					else:
						f_w.write("no;")
				else:
					cluster_dir = os.path.dirname(NoNewLine)
					if(os.path.isfile(cluster_dir + "/outgroup_seq_concat_shortdesc.fasta")):
						f_w.write(cluster_dir + "/outgroup_seq_concat_shortdesc.fasta;")
						if outgroup_seq_file==0:
							Nex_outgroupPath=(cluster_dir + "/outgroup_seq_concat_shortdesc.fasta")
							outgroup_seq_file=1
					else:
						f_w.write("no;")
				f_w.write(cluster_dir + "/mb_config1.nex;")
				f_w.write(cluster_dir + "/mb_final_seq1.nex;")
				f_w.write(cluster_dir + "/mb1.out;")
				f_w.write("no" +"\n")
		if outgroup_seq_file == 0:
			Nex_outgroupPath = "no"


	#default values:
	#nchains_val = 4
	#sampleFreq_val = 2000
	#burninFrac_val = 0.25
	#Checkfreq =  10000
	#10. nchains_val
	#11. samplefreq_val
	#12. burninFrac_val
	#13. checkFreq_val

	#Check if node dating is perfromed or not:
	# If so, add the following block:

	#constraint split_node = Areca_rheophytica Areca_concinna;
	#calibrate split_node = uniform(0.5,0.7);
	#prset topologypr = constraints(split_node);
	#prset nodeagepr = calibrated;


	if split_flag == 'no':
		if(DEBUG_PRINT == 1): print(("perl " + scripts_dir + "create_nexus_for_concat.pl " + concatDir + "cluster_to_concat_list.txt " +
		concatDir+genus_name+"-concat-aligned.fasta " + concatDir+genus_name+"-concat-aligned.report " + Nex_outgroupPath +
		""" "ngen=""" + Ngen_param[0] + """" """ + concatDir+"mb_config.nex " + concatDir+"mb_final_seq.nex " + concatDir+"mb.out " + "no " +
						str(Ngen_param[1]) + " " + str(Ngen_param[2]) + " " + str(Ngen_param[3]) + " " + str(Ngen_param[4]) +
						" " + clock_model + " " + relaxed_branch ))

		f_w.close()
		os.system("perl " + scripts_dir + "create_nexus_for_concat.pl " + concatDir + "cluster_to_concat_list.txt " +
		concatDir+genus_name+"-concat-aligned.fasta " + concatDir+genus_name+"-concat-aligned.report " + Nex_outgroupPath +
		""" "ngen=""" + Ngen_param[0] + """" """ + concatDir+"mb_config.nex " + concatDir+"mb_final_seq.nex " + concatDir+"mb.out " + "no " +
		str(Ngen_param[1]) + " " + str(Ngen_param[2]) + " " + str(Ngen_param[3]) + " " + str(Ngen_param[4]) +
		" " + clock_model + " " + relaxed_branch)
	else:
		print(Ngen_param[5])
		if(DEBUG_PRINT == 1): print(("perl " + scripts_dir + "create_nexus_for_concat.pl " + concatDir + "cluster_to_concat_list.txt " +
		concatDir+genus_name+"-concat-aligned.fasta " + concatDir+genus_name+"-concat-aligned.report " + Nex_outgroupPath +
		""" "ngen=""" + Ngen_param[0] + """" """ + concatDir+"mb_config.nex " + concatDir+"mb_final_seq.nex " + concatDir+"mb.out " + "no " +
						str(Ngen_param[1]) + " " + str(Ngen_param[2]) + " " + str(Ngen_param[3]) + " " + str(Ngen_param[4]) +
						" " + clock_model + " " + relaxed_branch + " " + Ngen_param[5] + "'"))

		f_w.close()
		os.system("perl " + scripts_dir + "create_nexus_for_concat.pl " + concatDir + "cluster_to_concat_list.txt " +
		concatDir+genus_name+"-concat-aligned.fasta " + concatDir+genus_name+"-concat-aligned.report " + Nex_outgroupPath +
		""" "ngen=""" + Ngen_param[0] + """" """ + concatDir+"mb_config.nex " + concatDir+"mb_final_seq.nex " + concatDir+"mb.out " + "no " +
		str(Ngen_param[1]) + " " + str(Ngen_param[2]) + " " + str(Ngen_param[3]) + " " + str(Ngen_param[4]) +
		" " + clock_model + " " + relaxed_branch + " " + Ngen_param[5] + "'")


	#Add Node dating test in case needed:
	if split_flag != 'no':
		for line in fileinput.FileInput(concatDir+"mb_config.nex",inplace=1):
			if "mcmc nruns" in line:
				line=line.replace(line,Ngen_param[5]+line)
				print(line)
			else:
				print(line.rstrip())

	#sys.exit()



	return

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def part1_create_jmt(working_dir,scripts_dir,fastaFilesList,standAlone_Flag):
# Create phylip files from fasta files for each concat cluster, according to the list at MSA output directory :
# GENERA/concat/fasta-files-to-concat.txt
#fastaFilesList = (MSA_OUTPUT_PATH + genus_name + "/concat/fasta-files-to-concat.txt")       # build path for fasta files to concat

	concat_dir = os.path.dirname(fastaFilesList)
	jobs_names_list=[]
	jobs_dirs_list=[]
	f_jobs = open(concat_dir + '/jmt_JobsList.txt','w')
	with open(fastaFilesList, "r") as pathFile:
		for FastaPath in pathFile:
			file_no_extension, file_extension = os.path.splitext(FastaPath)	# remove file extension
			phy_file = file_no_extension+PHY_EXT							# full path + file
			fasta_file=FastaPath.replace('\n','')							# Remove \n from path
			msa_phy_file=(os.path.basename(FastaPath)).replace('.fasta\n',PHY_EXT)			# fasta msa file name only !!!
			phy_file=file_no_extension+PHY_EXT
			# 1. Convert cluster fasta file to phylip:
			convert_fasta_to_phylip(fasta_file, phy_file)
			# Cluster dir path - folder only
			cluster_dir = os.path.dirname(FastaPath)
			jobs_dirs_list.append(cluster_dir)
			# 2. Create jModelTest - run job file
			jmt_file_path=file_no_extension + ".jmt"			# jmt output file path + file name
			jmt_JOB_file_path=file_no_extension + "_jmtJOB.sh"
			sh_filePath=jmt_JOB_file_path
			if(DEBUG_PRINT == 1): print (sh_filePath)
			#Create Genus Job names -
			cluster_dir_notFull = cluster_dir.rsplit('/',1)[1]
			working_dir_only = working_dir.rsplit('/',2)[1]
			job_name = '-'.join(["jmt", working_dir_only,cluster_dir_notFull])
			#Cancel qsub for standAlone version:
			jmodelTest_path = ploidb_config['diff_soft']['jmodelTest']
			if standAlone_Flag == 'off': #default
				job_filename = createJobFile.create_job_file(job_name, " ".join(["python " + scripts_dir +
								"run_jModelTest.py -p", cluster_dir+"/"+msa_phy_file, "-m", msa_phy_file,
								"-jm", jmodelTest_path, "-o", cluster_dir]), sh_filePath, cluster_dir)
				if(DEBUG_PRINT == 1): print("Job file name is : %s" % job_filename)
				os.system('qsub ' + job_filename)
				jobs_names_list.append(job_name)
				f_jobs.write(job_name+'\n')
			else:
				jmodelTest_path = ploidb_config['diff_soft']['jmodelTest']
				cmd_jmt = "python " + scripts_dir + "run_jModelTest.py -p " + cluster_dir + "/" + msa_phy_file + \
								" -m " + msa_phy_file + " -jm " + jmodelTest_path + " -o " +  cluster_dir
				os.system(cmd_jmt)

	print("JOB NAMES LIST *********************************************")
	print(jobs_names_list)
	f_jobs.close()
	return jobs_dirs_list

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def mb_UserModel_part12(working_dir,scripts_dir,fastaFilesList,User_Model):
# Create phylip files from fasta files for each concat cluster, according to the list at MSA output directory :
# Then create BLK files for MB run

	concat_dir = os.path.dirname(fastaFilesList)
	#jobs_names_list=[]
	#jobs_dirs_list=[]
	#f_jobs = open(concat_dir + '/jmt_JobsList.txt','w')
	with open(fastaFilesList, "r") as pathFile:
		for FastaPath in pathFile:
			file_no_extension, file_extension = os.path.splitext(FastaPath)	# remove file extension
			phy_file = file_no_extension+PHY_EXT							# full path + file
			fasta_file=FastaPath.replace('\n','')							# Remove \n from path
			msa_phy_file=(os.path.basename(FastaPath)).replace('.fasta\n',PHY_EXT)			# fasta msa file name only !!!
			phy_file=file_no_extension+PHY_EXT
			# 1. Convert cluster fasta file to phylip:
			convert_fasta_to_phylip(fasta_file, phy_file)
			cluster_dir= os.path.dirname(FastaPath)
			fixed_model_and_create_blk(phy_file,cluster_dir,User_Model)
			logger.debug("Create Blk file at: %s" %cluster_dir)
		logger.debug("Create Blk files for all clusters")
	return

#-----------------------------------------------------------------------------------------------------------------------

def Create_partitionFile(report_file,xml_partition_file):

	gene_start_end_dict={}
	with open(report_file,'r') as f_report:
		firstline = True
		for line in f_report:
			if firstline == True:
				firstline=False
				continue
			columns = line.split('\t')
			if (int(columns[1]) + 1) not in gene_start_end_dict.keys():
				gene_start_end_dict[int(columns[1]) + 1]=int(columns[2])
	f_partition = open(xml_partition_file,'w')
	gen_ind=1
	for key in gene_start_end_dict.keys():
		f_partition.write("DNA, gene_%d = %d-%d\n" %(gen_ind,key,gene_start_end_dict[key])) # DNA, gene_18S.mafft = 1-1836
		gen_ind+=1
	f_partition.close()

