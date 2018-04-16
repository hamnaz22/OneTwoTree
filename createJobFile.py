
__author__ = 'Shiran'

import argparse


def create_job_file(job_name, command, file_name, error_files_path):

	with open(file_name, "w") as handle:
		handle.write("#!/bin/tcsh\n\n")
		handle.write("#$ -N " + job_name + "\n")
		handle.write("#$ -S /bin/tcsh\n")
		handle.write("#$ -cwd\n")
		#handle.write("#$ -l itaym\n")
		handle.write("#$ -e " + error_files_path + "$JOB_NAME.$JOB_ID.ER\n")
		handle.write("#$ -o " + error_files_path + "$JOB_NAME.$JOB_ID.OU\n")
		#handle.write("module load python/python-3.3.0\n")
		#handle.write("module load perl/perl-5.16.3\n")
		handle.write("module load mrbayes/mrbayes_3.2.2\n")
		handle.write(command + "\n")
	return file_name

def create_job_file_TwoCmd(job_name, command1, command2, file_name, error_files_path):

	with open(file_name, "w") as handle:
		handle.write("#!/bin/tcsh\n\n")
		handle.write("#$ -N " + job_name + "\n")
		handle.write("#$ -S /bin/tcsh\n")
		handle.write("#$ -cwd\n")
		#handle.write("#$ -l itaym\n")
		handle.write("#$ -e " + error_files_path + "$JOB_NAME.$JOB_ID.ER\n")
		handle.write("#$ -o " + error_files_path + "$JOB_NAME.$JOB_ID.OU\n")
		#handle.write("module load python/python-3.3.0\n")
		#handle.write("module load perl/perl-5.16.3\n")
		handle.write("module load mrbayes/mrbayes_3.2.2\n")
		handle.write(command1 + "\n")
		handle.write(command2 + "\n")
	return file_name

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Create jobs for queue')
	parser.add_argument('--shFile','-s', help='sh filename', required=True)
	parser.add_argument('--command','-c', help='Command to execute', required=True)
	parser.add_argument('--jobName','-j', help='Job name', required=True)

	args = parser.parse_args()
	create_job_file(args.jobName, args.command, args.shFile)
