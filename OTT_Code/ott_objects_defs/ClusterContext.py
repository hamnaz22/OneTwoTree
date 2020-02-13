import os
import re
import operator
from ott_objects_defs.OutGroupSelection import OutGroupSelection
from ott_objects_defs.handleMultAccessions import get_seq_with_max_average_blast_score
from ott_objects_defs.ploidbCommon import *
from Bio import SeqIO


__author__ = 'moshe'

"""
Class for storing a single cluster data
"""

DEBUG_FLAG = 0

class SeqType(object):
	Nuc, Mt, Chloroplast = range(3)


class ClusterContext(object):
	def __init__(self, index, all_seqs_fasta_filename, genus_cluster):
		"""
		@type genus_cluster: PloiDbContext
		"""
		self.id = genus_cluster.id
		self.genus_working_dir = genus_cluster.working_dir
		self.cluster_id = ""
		self.cluter_work_dir = ""
		self.all_seqs_fasta_filename = ""
		self.outgroup_seq_filename = ""
		self.all_seqs_with_outgroup_filename = ""
		self.is_its_cluster = False
		self.locus_tree_outgroup_selection = OutGroupSelection()
		self.cluster_desc = ""
		self.cluster_type = None
		self.is_mult_acc_init = False
		self.is_db_init = False
		self.index = index
		self.is_used_for_locus_tree = False
		self.is_used_for_concat_tree = False
		self.number_of_different_species = None
		self.original_fasta_filename = all_seqs_fasta_filename
		#self.gb_seqs_dict = None
		self.data_matrix_size = None
		self.seq_length = None
		self.estimated_seq_length = None
		self.no_of_outgroup_sequence_for_concat = -1
		self.cluster_seq_to_blast = None
		self.gene_name = None
		self.rep_blast_seq_id = None

		self.cluster_taxid_accession_dict = {}

		self.cluster_type_dict = {}
		#self.gb_seqs_dict = genus_cluster.get_#()
		self.init_dirs(genus_cluster)

	#def __getstate__(self):
	#	odict = self.__dict__.copy()
	#	if 'gb_seqs_dict' in odict:
	#		del odict['gb_seqs_dict']
	#	return odict
#
	#def get_gb_seqs_dict(self):
	#	if self.gb_seqs_dict is None:
	#		self.gb_seqs_dict = get_indexed_seq_hash(ott_config['general']['genbankDir'])
#
	#	return self.gb_seqs_dict

	def init_dirs(self, genus_cluster):
		self.cluter_work_dir = os.path.join(genus_cluster.working_dir, str(self.index))
		self.cluster_concat_work_dir = os.path.join(genus_cluster.concat_workdir, str(self.index))
		self.all_seqs_fasta_filename = os.path.join(self.cluter_work_dir, "original-seqs.fasta")
		self.all_seqs_for_blast_fasta_filename = os.path.join(self.cluter_work_dir, "original-seqs-for-blast.fasta")
		self.init_mb_dirs(genus_cluster.mrbayes_input_dir, genus_cluster.mrbayes_output_dir)
		self.rewrite_names_dict_filename = os.path.join(self.cluter_work_dir, "rewritten-names.json")

		self.blast_results_filename = os.path.join(self.cluter_work_dir, "blast_all-v-all_total.blastn")
		self.multiple_accessions_work_dir = os.path.join(self.cluter_work_dir, "multi-acc-temp")
		create_dir_if_not_exists(self.multiple_accessions_work_dir)

		self.outgroup_seq_full_filename = self.cluter_work_dir + "/outseq-fullname.fasta"
		self.outgroup_seq_filename = self.cluter_work_dir + "/outseq-organism.fasta"
		self.all_seqs_fasta_filename_no_multiple_accessions_temp = self.cluter_work_dir + "/seqs-no-multiple-accessions-temp.fasta"
		self.all_seqs_fasta_filename_no_multiple_accessions = self.cluter_work_dir + "/seqs-no-multiple-accessions.fasta"
		self.all_seqs_with_outgroup_filename = self.cluter_work_dir + "/seqs-with-out-group.fasta"
		self.all_seqs_short_desc = self.cluter_work_dir + "/all-seqs-short-desc.fasta"
		self.reverse_seqs_in_fasta = self.cluter_work_dir + "/seqs-after-rev.fasta"
		self.all_seqs_aligned_fasta_filename = self.cluter_work_dir + "/seqs-aligned.fasta"
		self.all_seqs_for_tree_creation_fasta_filename = self.cluter_work_dir + "/seqs-organism.fasta"
		self.raxml_work_dir = os.path.join(self.cluter_work_dir, "raxml")


		# Concat files
		self.outgroup_sequence_for_concat_filename = self.cluster_concat_work_dir + "/outgroup_seq_concat.fasta"
		self.outgroup_sequence_for_concat_short_desc_filename = self.cluster_concat_work_dir + "/outgroup_seq_concat_shortdesc.fasta"
		self.all_seqs_with_concat_outgroup_filename = self.cluster_concat_work_dir + "/seqs-with-out-group-concat.fasta"
		self.all_seqs_concat_short_desc = self.cluster_concat_work_dir + "/all-seqs-cocnat-short-desc.fasta"
		self.reverse_seqs_cocnat_in_fasta = self.cluster_concat_work_dir + "/seqs-after-rev-concat.fasta"
		self.all_seqs_concat_aligned_fasta_filename = self.cluster_concat_work_dir + "/seqs-aligned-concat.fasta"
		self.all_seqs_for_concat_tree_creation_fasta_filename = self.cluster_concat_work_dir + "/seqs-organism-concat.fasta"
		self.raxml_concat_work_dir = os.path.join(self.cluster_concat_work_dir, "raxml")

		# Guidance
		self.guidance_work_dir1 = os.path.join(self.cluter_work_dir, "guidance1")
		self.guidance_work_dir2 = os.path.join(self.cluter_work_dir, "guidance2")
		self.guidance_concat_work_dir1 = os.path.join(self.cluster_concat_work_dir, "guidance1")
		self.guidance_concat_work_dir2 = os.path.join(self.cluster_concat_work_dir, "guidance2")

		# NOG = no out group
		self.all_seqs_concat_aligned_nog_fasta_filename = self.cluster_concat_work_dir + "/seqs-aligned-concat-nog.fasta"
		self.all_seqs_for_concat_tree_creation_nog_fasta_filename = self.cluster_concat_work_dir + "/seq-organism-concat-nog.fasta"

	def init_cluster(self, genus_context):
		if not os.path.exists(self.cluter_work_dir):
			logger.debug("Creating cluster work dir %s" % self.cluter_work_dir)
			os.makedirs(self.cluter_work_dir)
		if not os.path.exists(self.cluster_concat_work_dir):
			logger.debug("Creating cluster work dir %s" % self.cluster_concat_work_dir)
			os.makedirs(self.cluster_concat_work_dir)

		logger.debug("Writing sorted %s to %s - as the orginal cluster seqs file" % (
			self.original_fasta_filename, self.cluter_work_dir))
		#Check for ITS:
		#if 'oneSeqPerSpecies' in self.original_fasta_filename:
		#	with open(self.cluter_work_dir + '/ITS_CLUSTER', 'w') as f_flag_its:
		#		f_flag_its.close()
		all_seqs_orig = list(SeqIO.parse(self.original_fasta_filename, "fasta"))
		self.num_of_cluster_seqs = len(all_seqs_orig)

		all_seq_sorted = sorted(all_seqs_orig, key=lambda next_seq: next_seq.id)
		with open(self.all_seqs_fasta_filename, "w") as handle:
			SeqIO.write(all_seq_sorted, handle, "fasta")

		#Cast of ITS cluster
		logger.debug("INSIDE init_cluster\n")
		logger.debug("self.original_fasta_filename : %s\n" %self.original_fasta_filename)
		logger.debug("self.cluter_work_dir : %s\n" %self.cluter_work_dir)
		logger.debug("self.cluster_concat_work_dir : %s\n" %self.cluster_concat_work_dir)
		if 'oneSeqPerSpecies' not in self.original_fasta_filename:
			self.init_blast(genus_context)
			self.init_cluster_id()
		else:
			ITS_options = ['ITS_cluster_2','ITS_cluster_3','ITS_cluster_4','ITS_cluster_5']
			#/groups/itay_mayrose/michaldrori/OTT_TestingDir/ITS_clustering/cluster_its
			Other_ITS_fount=0
			for its_dir in ITS_options:
				if its_dir in self.original_fasta_filename:
					Other_ITS_fount = 1
					logger.debug("ITS found more clusters: %s" % (self.cluter_work_dir + '/' + its_dir))
					# get ITS number from path and add to ID:
					r = re.compile('ITS_cluster_(\d)')
					m = r.search(its_dir)
					if m:
						logger.debug("ITS Num cluster %s" % m.group(1))
						ITS_num = str(m.group(1))
						with open(self.cluter_work_dir + '/ITS_CLUSTER_' + ITS_num, 'w') as f_flag_its:
							f_flag_its.close()
						self.init_blast(genus_context)
						self.cluster_desc = 'ITS' + ITS_num
						self.cluster_id = 'ITS' + ITS_num
			if 'cluster_its' in self.original_fasta_filename and Other_ITS_fount == 0:
				logger.debug("ITS cluster main: %s" % self.cluter_work_dir)
				with open(self.cluter_work_dir + '/ITS_CLUSTER', 'w') as f_flag_its:
					f_flag_its.close()
				self.init_blast(genus_context)
				self.cluster_desc = 'ITS'
				self.cluster_id = 'ITS'

		self.init_cluster_type()


	def init_blast(self, genus_context):
		exec_blast_all_vs_all(self.all_seqs_fasta_filename, self.all_seqs_for_blast_fasta_filename,
							  self.blast_results_filename, genus_context.id)


	def get_cluster_representative(self):
		allRecords = list(SeqIO.parse(self.all_seqs_fasta_filename, "fasta"))
		self.num_of_cluster_seqs = len(allRecords)
		return allRecords[0]


	def get_cluster_seq_to_blast(self):
		cluster_seq_to_blast = self.cluster_seq_to_blast
		if cluster_seq_to_blast is None:
			cluster_seq_to_blast = get_seq_with_max_average_blast_score(self.all_seqs_fasta_filename,
																		self.blast_results_filename)
			#seqid = getPropertyFromFastaSeqHeader(cluster_seq_to_blast.description, "seqid")
			seqid = getPropertyFromFastaSeqHeader(cluster_seq_to_blast.description, "gi")
			logger.debug("Representative seq for cluster BLAST is %s" % seqid)
			#self.rep_blast_seq_id = seqid
		#else:
			#self.rep_blast_seq_id = cluster_seq_to_blast
		return cluster_seq_to_blast


	def init_cluster_type(self):
		self.cluster_type = SeqType.Nuc
		self.gene_name = None

		gene_names_to_counts = dict()
		mt_count = 0
		chloroplast_count = 0
		other_count = 0
		none_count = 0


		#gb_seqs_dict = self.get_gb_seqs_dict()
		all_seqs = list(SeqIO.parse(self.all_seqs_fasta_filename, "fasta"))

		dir_path_ = os.path.dirname(os.path.abspath(self.all_seqs_fasta_filename))

		genes_list=[]
		for seq in all_seqs:
			one_gene = ''
			no_brack=''
			#seq_id = getPropertyFromFastaSeqHeader(seq.description, "seqid")  #try to replace with gi sionce the seqid sometimes is chopped of because of headline length
			seq_id = getPropertyFromFastaSeqHeader(seq.description, "gi")
			##if seq_id not in ##:
			##	#logger.warning("seqid %s was not found in dict" % seq_id)
			##	logger.warning("seqid %s was not found in dict" % seq_id)
			##	return
			##seq_gb = gb_seqs_dict[seq_id]
			logger.debug("Checking for key cluster type in seqid-%s desc-%s" % (seq_id,seq.description))
			if "chloroplast" in seq.description:
				chloroplast_count += 1
			elif "mitoch" in seq.description:
				mt_count += 1
			else:
				other_count += 1

			#logger.debug(" @@@ SeqID equals to rep_blast_seq_id: %s" %self.rep_blast_seq_id)
			if seq_id == self.cluster_id:
				#logger.debug(" @@@ SeqID equals to rep_blast_seq_id: %s" %self.rep_blast_seq_id)

				if "chloroplast" in seq.description:
					chloroplast_count += 1
				elif "mitoch" in seq.description:
					mt_count += 1
				else:
					other_count += 1

				##for feature in seq.features:
				##	if DEBUG_FLAG == '1': logger.debug('init_cluster_type - print feature typefor debug:')
				##	if DEBUG_FLAG == '1': logger.debug(feature)
				##	if feature.type == "source":
				##		organelles = feature.qualifiers.get("organelle")
				##		if organelles is None:
				##			none_count += 1
				##		else:
				##			for organelle in organelles:
				##				if DEBUG_FLAG == '1': logger.debug("organelle is %s" % organelle)
				##				if DEBUG_FLAG == '1': logger.debug("MDebug: organelle is %s" %organelle)
				##				if organelle.find("chloroplast") != -1 or organelle.find("plastid") != -1:   # Added PlastID to be recognized as chloroplast as well !!! (Michal D.)
				##					chloroplast_count += 1
				##				elif organelle.find("mitoch") != -1 or organelle.find("mt") != -1:
				##					mt_count += 1
				##				if DEBUG_FLAG == '1': logger.debug("Found organelle %s" % organelle)
				##				else:
				##					other_count += 1
				##	elif feature.type == "gene":
				##		one_gene = feature.qualifiers.get("gene")
				##		if one_gene:
				##			for item in one_gene:
				##				no_brack = re.sub(r' \(.*?\)', '', item)
				##				#if no_brack not in genes_list:
				##				found_CaseIdent=0
				##				for gene_in_list in genes_list:
				##					if gene_in_list.upper() == no_brack.upper():
				##						found_CaseIdent=1
				##				if found_CaseIdent != 1:
				##					genes_list.append(no_brack)
				##			#if no_brack not in genes_list:
				##			#	genes_list.append(no_brack)
				##		#genes_str = "|".join(genes)
##
				##		#if genes is not None and len(genes) > 0:
				##		#	#logger.debug("Number of gene qualifiers is %d" % len(genes))
				##		#	genes = sorted(genes)
				##		#	genes_str = "|".join(genes)
				##		#	#logger.debug("Found the following genes %s" % genes_str)
				##		#	if not genes_str in gene_names_to_counts:
				##		#		gene_names_to_counts[genes_str] = 0
				##		#	gene_names_to_counts[genes_str] += 1
				##	else:
				##		pass

					#if genes_str != '':
					#	logger.debug("Found the following genes %s (seqid %s)" % (genes_str,seq_id))
					#	if not genes_str in gene_names_to_counts:
					#		gene_names_to_counts[genes_str] = 0
					#	gene_names_to_counts[genes_str] += 1

					#logger.debug("Ignoring feature %s" % feature.type)

		genes_list = list(set(genes_list))

		if genes_list:
			self.gene_name = '|'.join(genes_list)
		logger.info("none_count=%d chloroplast_count=%d mt_count=%d other_count=%d" % (
			none_count, chloroplast_count, mt_count, other_count))
		total_seqs = len(all_seqs)
		# Set Type of Clsuter according to majority rule:
		nuc_count = none_count + other_count

		ClustType_str = "NotDefined"
		if chloroplast_count == mt_count == nuc_count:
			ClustType_str = "NotDefined"
		else:
			#Check which one has the maximum value:
			max_type = max(chloroplast_count, mt_count, nuc_count)
		if max_type == chloroplast_count:
			logger.debug("seqType Chloroplast: %d out of %d " %(chloroplast_count,total_seqs))
			self.cluster_type = SeqType.Chloroplast
			ClustType_str = "cpDNA"
		#elif 1.0 * mt_count / total_seqs > 0.9:
		elif max_type == mt_count:
			logger.debug("seqType Mt: %d out of %d " %(mt_count,total_seqs))
			self.cluster_type = SeqType.Mt
			ClustType_str = "mtDNA"
		else:
			logger.debug("seqType Nuc: %d out of %d "%(nuc_count,total_seqs))
			self.cluster_type = SeqType.Nuc
			ClustType_str = "NUC"
		#Create cluster Type file and update context dict:
		Cluster_indx = os.path.basename(os.path.normpath(dir_path_))
		self.cluster_type_dict[Cluster_indx] = ClustType_str
		f_temp = open(dir_path_ + "/ClustType_" + ClustType_str,'w')
		f_temp.close()

	def get_cluster_identifier(self):
		desc = self.gene_name if self.gene_name is not None else self.cluster_desc
		return "%s-%s" % (str(self.index), desc)


	def init_cluster_id(self):
		cluster_rep = self.get_cluster_seq_to_blast()
		self.cluster_id = getPropertyFromFastaSeqHeader(cluster_rep.description, "gi")
		logger.debug(" Inside init_cluster_id: cluster rep is %s" % self.cluster_id)
		self.cluster_desc = getPropertyFromFastaSeqHeader(cluster_rep.description, "description")
		#check if descriptuion is None
		if not self.cluster_desc:
			self.cluster_desc = 'None'


	def init_mb_dirs(self, mrbayes_input_dir_for_all_genus, mrbayes_output_dir_for_all_genus):
		self.mrbayes_input_dir = mrbayes_input_dir_for_all_genus + "/" + str(self.index)
		self.mrbayes_output_dir = mrbayes_output_dir_for_all_genus + "/" + str(self.index)

		# Initalizing the paths for mr AIC / mr Bayes output for the trees per locus
		self.mb_config_file = os.path.join(self.mrbayes_input_dir, "mb_config1.nex")
		self.mb_config_file1 = os.path.join(self.mrbayes_input_dir, "mb_config1.nex1")
		self.mb_config_file2 = os.path.join(self.mrbayes_input_dir, "mb_config1.nex2")
		self.mb_seq_file = os.path.join(self.mrbayes_input_dir, "mb_final_seq1.nex")
		self.mb_out_file = os.path.join(self.mrbayes_output_dir, "mb1.out")
		self.mb_out_file1 = os.path.join(self.mrbayes_output_dir, "mb1.out1")
		self.mb_out_file2 = os.path.join(self.mrbayes_output_dir, "mb1.out2")

		self.mb_out_consensus_file = os.path.join(self.mrbayes_output_dir, "mb1.out.con.tre")
		self.mb_out_consensus_file_newick = os.path.join(self.mrbayes_output_dir, "mb1.out.con.newick")


		# Initalizing the paths for mr AIC / mr Bayes output for the concat trees
		self.mb_concat_config_file = os.path.join(self.cluster_concat_work_dir, "mb_config_concat.nex")
		self.mb_concat_seq_file = os.path.join(self.cluster_concat_work_dir, "mb_final_seq_concat.nex")
		self.mb_concat_out_file = os.path.join(self.cluster_concat_work_dir, "mb-concat.out")

		# Just for cache:
		self.mraic_ouput_file = os.path.join(self.mrbayes_input_dir, "seqs-organism-concat_phy.phy.MrAIC.txt")


	def get_newick_concat_consensus_tree(self):
		return get_newick_tree_str(self.mb_out_consensus_file_newick)


	def get_view_name(self):
		return "c" + (self.cluster_id).replace('.','_')


	def get_list_of_species(self):
		return get_list_of_species(self.all_seqs_fasta_filename)


	def get_number_of_different_species(self):
		if self.number_of_different_species is None:
			logger.debug("Getting number_of_different_species from %s" % self.all_seqs_fasta_filename)
			#self.number_of_different_species = ploidbCommon.get_number_of_different_species(self.all_seqs_fasta_filename)
			self.number_of_different_species = get_number_of_different_species(self.all_seqs_fasta_filename)
		return self.number_of_different_species


	def get_cluster_length(self, estimated=False):
		seq_length = self.seq_length
		if seq_length is None:
			#If ITS_CLUSTER calc the median length of combined ?????????
			if os.path.exists(self.all_seqs_aligned_fasta_filename):
				all_seqs = list(SeqIO.parse(self.all_seqs_aligned_fasta_filename, "fasta"))
				self.seq_length = len(all_seqs[0])
				seq_length = self.seq_length
			else:
				if estimated and os.path.exists(self.all_seqs_fasta_filename):
					all_seqs = list(SeqIO.parse(self.all_seqs_fasta_filename, "fasta"))
					all_seq_sorted = sorted(all_seqs, key=lambda next_seq: len(next_seq))
					median_index = int(len(all_seq_sorted) / 2)
					#seq_length = len(all_seqs[median_index]) # Bug !!! fixed March17
					seq_length = len(all_seq_sorted[median_index])
					self.estimated_seq_length = seq_length
				else:
					raise Exception(
						"Couldn't find any file to determine the cluster length. Setting length = 0 for cluster %s (index=%s)" % (
							self.cluster_id, str(self.index)))
		return seq_length


	def ITS_blast_all(self,ploidb_context):
		#its_cluster_dir
		all_seqs_its = ploidb_context.cluter_its_dir + '/oneITSTypePerSpecies.fasta'
		if not os.path.exists(self.blast_results_filename):
			exec_blast_all_vs_all(all_seqs_its, self.all_seqs_for_blast_fasta_filename,
							  self.blast_results_filename, ploidb_context.id)
		return

	#This ratio is used for cluster selection process:
	# L_ratio = number of data characters in alignment decided by the total length of the final MSA:
	def get_cluster_L_ratio(self,ploidb_context):
		# IF ITS then we need to run blast_all_vs_all:
		ITS_options = ['ITS_CLUSTER','ITS_CLUSTER_2', 'ITS_CLUSTER_3', 'ITS_CLUSTER_4', 'ITS_CLUSTER_5']
		for its_op in ITS_options:
			if os.path.exists(self.cluter_work_dir + '/' + its_op):
				self.ITS_blast_all(ploidb_context)

		#Since we don't have the msa as this point in the pipeline we need some estimation of the msa quality.
		#For that we'll use the blast_all_v_all we run for each cluster (for handling the multiple accessions):

		seq_id_length_dict = dict()
		seq_id_list = []
		if os.path.exists(self.blast_results_filename):
			with open(self.blast_results_filename) as f:

				dr = csv.DictReader(f, delimiter='\t',
									fieldnames=['query_id', 'subject_id', 'pct_identity', 'align_len', 'mismatches',
												'gap_openings', 'q_start', 'q_end', 's_start', 's_end', 'eval',
												'bitscore'])

				logger.info("Set L_ratio for msa quality param for cluster selection:")
				list_pct_identity = []
				dict_pct_identity = dict()
				dict_qid_pct_identity = dict()
				for next_row in dr:
					qid = next_row['query_id']
					sid = next_row['subject_id']
					#Collect sequences length and calc total align_len (for cases where the alignment is cut
					# into more than one section)
					if qid == sid:
						if qid in seq_id_length_dict.keys():
							seq_id_length_dict[qid] += float(next_row['align_len'])
						else:
							seq_id_length_dict[qid] = float(next_row['align_len'])
						#logger.debug("seq_id_length_dict[qid] = %f" % (seq_id_length_dict[qid]))

			with open(self.blast_results_filename) as f2:
				dr2 = csv.DictReader(f2, delimiter='\t',
								fieldnames=['query_id', 'subject_id', 'pct_identity', 'align_len', 'mismatches',
											'gap_openings', 'q_start', 'q_end', 's_start', 's_end', 'eval',
											'bitscore'])

				for next_row in dr2:
					ref_length=0
					qid = next_row['query_id']
					sid = next_row['subject_id']

					if qid != sid:
						qid_sid_pair = qid + '_' + sid
						if qid not in seq_id_list:
							seq_id_list.append(qid)
						ref_length = seq_id_length_dict[qid]
						# 357/940*0.83 +
						prec_align = float(next_row['align_len'])/seq_id_length_dict[qid]*float(next_row['pct_identity'])
						if qid_sid_pair in dict_pct_identity.keys():
							dict_pct_identity[qid_sid_pair] += prec_align
						else:
							dict_pct_identity[qid_sid_pair] = prec_align


				for qid in seq_id_list:
					for key in dict_pct_identity.keys():
						if qid in key:
							if qid in dict_qid_pct_identity.keys():
								dict_qid_pct_identity[qid].append(dict_pct_identity[key])
							else:
								dict_qid_pct_identity[qid] = [(dict_pct_identity[key])]


				for qid_seq in dict_qid_pct_identity.keys():
					prsct_ident = sum(dict_qid_pct_identity[qid_seq])/len(dict_qid_pct_identity[qid_seq])
					logger.debug("prsct_ident of seq id %s is %f\n" % (qid_seq,prsct_ident))
					list_pct_identity.append(prsct_ident)

			L_ratio = sum(list_pct_identity)/len(list_pct_identity)
			logger.info("L_ratio is %f" %L_ratio)
			return L_ratio
		else:
			raise Exception("Couldn't find blast_all_v_all for cluster %s (index=%s)" % (
					self.cluster_id, str(self.index)))



	def get_data_matrix_size(self, estimated=False):
		data_matrix_size = self.data_matrix_size
		if data_matrix_size is None:
			# cluster_length = self.get_cluster_length(estimated)
			#data_matrix_size = ((self.get_number_of_different_species())**2) * cluster_length
			data_matrix_size = self.get_number_of_different_species()
			if not estimated:
				self.data_matrix_size = data_matrix_size

		return data_matrix_size


	def __str__(self):
		#logger.warning("LookAtMe -> Cluster index %s, Cluster id %s" %(str(self.index), str(self.cluster_id)))
		#logger.warning("LookAtMe -> Cluster_desc %s" %(str(self.cluster_desc)))
		cluster_str = "[Cluster index=%s id=%s" % (str(self.index), self.cluster_id)
		try:
			if hasattr(self, "gene_name") and self.gene_name is not None:
				cluster_str += " gene=" + self.gene_name
			else:
				cluster_str += " desc=" + self.cluster_desc
		except AttributeError as e:
			logger.warning("AttributeError caught while writing cluster_str - probably genus_name doesn't exist")

		cluster_str += "]"
		return cluster_str
