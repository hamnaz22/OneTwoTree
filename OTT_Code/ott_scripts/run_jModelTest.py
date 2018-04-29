import os
import argparse
#from parse_jModelTest2_file import *


JMT_JAR_FILE = "jModelTest.jar"
SEP = "/"

DEBUG_PRINT = 1
PRE_TREE_DIR = 'PreTree_Dir'
PHY_EXT = '_phy.phy'
JMT_EXT = '_phy.phy.jmt'
FASTA_FORMAT = "fasta"
PHYLIP_FORMAT = "phylip-relaxed"


def run_jmt(msa_file_path, output_directory, jmodeltest_path,filename):
	output_file_path = create_output_file(SEP.join([output_directory, filename + ".jmt"]))
	if output_file_path == -1:
		return
	ch_dir = SEP.join([jmodeltest_path, JMT_JAR_FILE])
	try:
		os.system("java -jar " + ch_dir + " -d " + "\"" + msa_file_path + "\"" + \
				  " -g 4 -i -f -AIC -BIC -AICc -p -v -t BioNJ -o " + "\"" + output_file_path + "\"") \
		# " -tr 1 -g 4 -i -f -AIC -BIC -AICc -p -v -t BioNJ -o " + "\"" + output_file_path + "\"")	# -tr (trying to solve thred issue

		#Check if Jmodel file completed with models:: check for "::Best Models::"
		with open (output_file_path ,'r') as jmt_file:
			for line in jmt_file:
				if '::Best Models::' in line:
					jmt_dir_name = os.path.dirname(output_file_path)
					jmt_done_f = open(jmt_dir_name + '/JMT_DONE','w')		# Indication that current jmt is done
					jmt_done_f.close()
					return
	except Exception as exc:
		print (exc)
	return


def create_output_file(current_output_file_path):
	#don't overwrite a former output file
	#if os.path.exists(current_output_file_path):
	#	return -1
	fp = open(current_output_file_path, "w")
	fp.close()
	return current_output_file_path


if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Runs jMT2')
#	parser.add_argument('--genus_name', '-n', help='genus_name')
	parser.add_argument('--msa_file_name', '-m', help='msa_filename')
	parser.add_argument('--msa_file_path', '-p', help='msa_filepath')
	parser.add_argument('--jmodeltest_path', '-jm', help='JmodelTest_path')
	parser.add_argument('--output_dir', '-o', help='output_dir')
	args = parser.parse_args()

	# Run jModelTest:

	msa_file_path = args.msa_file_path
	print("CHECK THIS PATH: " + msa_file_path)
	if os.path.exists(msa_file_path):
		run_jmt(msa_file_path, args.output_dir, args.jmodeltest_path, args.msa_file_name)
