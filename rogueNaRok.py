import os
import re
from ploidbCommon import logger

__author__ = 'moshe'


# Rewrite MB config in the new location
# TODO: change this a little so that there is also a function for old conf file and new conf file (not based on dir)
def rewrite_mb_config(mb_config_filename, mb_new_out_dir,species_names_to_filter,new_mb_config_filename = None):
	mb_config_head, mb_config_tail = os.path.split(mb_config_filename)
	if new_mb_config_filename is None:
		new_mb_config_filename = os.path.join(mb_new_out_dir, mb_config_tail)
	mcmc_line_pattern = r"(\s*mcmc.*file=)(.*)/(.*)(;)"
	execute_line_pattern = r"(\s*execute\s*)(.*)/(.*)(;)"
	original_seq_nexus_file = None
	original_out_file = None
	with open(new_mb_config_filename, "w") as new_mb_conf_handle, open(mb_config_filename, "r") as old_mb_conf_handle:
		for line in old_mb_conf_handle:
			newline = line
			m = re.match(execute_line_pattern, line)
			if m:
				original_seq_nexus_file = m.group(2) + "/" + m.group(3)
				logger.info("Original nexus to exec is in %s" % original_seq_nexus_file)
				newline = re.sub(execute_line_pattern, r'\1%s/\3\4' % mb_new_out_dir, newline)
			m = re.match(mcmc_line_pattern, line)
			if m:
				original_out_file = m.group(2) + "/" + m.group(3)
				logger.info("Original out file is in %s" % original_out_file)
				newline = re.sub(mcmc_line_pattern, r'\1%s/\3\4' % mb_new_out_dir, newline)

			for s in species_names_to_filter:
				newline = newline.replace("%s " % s,"")
			new_mb_conf_handle.write(newline)

	return original_out_file, original_seq_nexus_file,new_mb_config_filename


def rewrite_mb_seq_file(new_mb_seq_filename, original_seq_nexus_file, species_names_to_filter):
	dim_line_pattern=r"(dimensions ntax)=(\d+)( nchar=)(\d+);"
	with open(new_mb_seq_filename, "w") as new_mb_seq_handle, open(original_seq_nexus_file,
																   "r") as original_mb_seq_handle:
		for line in original_mb_seq_handle:
			newline = line
			# Updating the number of species in the nexus file
			m = re.match(dim_line_pattern, newline)
			if m:
				old_num_of_species = float(m.group(2))
				new_num_of_species = old_num_of_species - len(species_names_to_filter)
				newline = re.sub(dim_line_pattern, r'\1=%d\3\4;' % new_num_of_species, newline)
				logger.info("Updating line %s to %s" % (line,newline))
				logger.info("Updated the number of species from %d to %d" % (old_num_of_species,new_num_of_species))

			# Removing the lines of the species_names_to_filter
			species_name_and_seq = newline.split("\t")
			matched_species = [s for s in species_names_to_filter if s == species_name_and_seq[0]]
			if len(matched_species) == 0:
				new_mb_seq_handle.write(newline)
			else:
				logger.info("Filtering %s from %s" % (",".join(matched_species), new_mb_seq_filename))

		#head, tail = os.path.split(new_mb_seq_filename)
		#(filename, ext) = os.path.splitext(tail)
		#before_align_fasta = os.path.join(head, filename + "-prealign" + ".fasta")
		#after_align_fasta = os.path.join(head, filename + "-postalign" + ".fasta")
		#
		#count = AlignIO.convert(new_mb_seq_filename_before_align, "nexus", before_align_fasta, "fasta", generic_dna)
		#logger.debug("Converted %i records from %s to %s" % (count,new_mb_seq_filename_before_align,before_align_fasta))
		#align_seq_file(before_align_fasta, after_align_fasta)
		#count = AlignIO.convert(after_align_fasta, "fasta", new_mb_seq_filename, "nexus", generic_dna)
		#logger.debug("Converted %i records from %s to %s" % (count,after_align_fasta,new_mb_seq_filename))

