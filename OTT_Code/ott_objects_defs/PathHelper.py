import fnmatch
import os
from ott_objects_defs.ClusterContext import SeqType

__author__ = 'Moshe'


class PathHelper:

	PICKLED_FILE_EXTENSION = ".pkl"

	def __init__(self, working_dir, id):
		self.id = id
		self.working_dir = working_dir + '/'
		self.cluster_method = "orthomclits"
		self.pickled_name = working_dir + "/" + self.id + "-" + self.cluster_method + self.PICKLED_FILE_EXTENSION
		self.pickled_name_after_trees = working_dir + "/" + self.id + "-" + "-aftertrees" + self.PICKLED_FILE_EXTENSION
		self.pickled_json_name = working_dir + "/" + self.id + "-" + self.cluster_method + ".json"
		self.pickled_json_name_after_trees = working_dir + "/" + self.id + "-aftertrees" + ".json"

		self.pickled_name_at_trees = working_dir + "/" + self.id + "-" + "-attrees" + self.PICKLED_FILE_EXTENSION
		self.pickled_json_name_at_trees = working_dir + "/" + self.id + "-attrees" + ".json"

		self.pre_tree_submitted =  os.path.join(working_dir, "submit-create-tree.sh")
		self.pre_tree_started =  os.path.join(working_dir, "pre_tree_started")
		self.pre_tree_stopped =  os.path.join(working_dir, "pre_tree_stopped")
		self.pre_tree_done =  os.path.join(working_dir, "pre_tree_done")
		self.create_tree_submitted =  os.path.join(working_dir, "submit-run-mrbayes.sh")
		self.create_tree_started =  os.path.join(working_dir, "create_tree_started")
		self.create_tree_done =  os.path.join(working_dir, "create_tree_done")
		self.db_file = os.path.join(working_dir,"ploidb.db")
		self.names_to_resolve = os.path.join(working_dir,"names_to_resolve.csv")

		self.synonym_names = os.path.join(working_dir,"synonym_names.csv")

		self.name_resolve_results = os.path.join(working_dir,"resolved_names.csv")
		self.ploidb_db_file = os.path.join(working_dir,"../../ploidb-full.db")
		self.fasta_seq_filename = os.path.join(working_dir,self.id + "-allseq" + ".fasta")
		self.fasta_seq_org_names_filename = os.path.join(working_dir,self.id + "-orgnames-allseq" + ".fasta")
		self.gb_seq_filename = os.path.join(working_dir,self.id + "-allseq" + ".gb")
		self.outgroup_workdir = os.path.join(working_dir,"outgroup-processing")
		self.mrbayes_input_dir = working_dir
		self.mrbayes_output_dir = working_dir
		self.clustering_dir = os.path.join(working_dir,"clustering")
		self.blast_results_filename = self.clustering_dir + "/blast_all-v-all_output.blastn"
		self.clustering_dir_inner = self.clustering_dir + "/inner"
		self.clustering_results_dir = os.path.join(working_dir,"cluster-results")
		self.cluter_its_dir = os.path.join(working_dir, "cluster_its")
		self.orthomcl_mcl_input = os.path.join(self.clustering_dir,"mclInput")


		#
		# This part is for the concat part
		#
		self.mrbayes_concat_work_dir = os.path.join(self.mrbayes_input_dir,"concat")
		self.raxml_concat_work_dir = os.path.join(self.mrbayes_concat_work_dir,"raxml")
		self.clusters_to_concat_filename = os.path.join(self.mrbayes_concat_work_dir,"cluster_to_concat_list.txt")
		self.concat_workdir = self.mrbayes_concat_work_dir
		self.concat_seqs_fasta_filename = os.path.join(self.mrbayes_concat_work_dir , self.id + "-concat-aligned.fasta")
		self.concat_seqs_report_filename = os.path.join(self.mrbayes_concat_work_dir , self.id + "-concat-aligned.report")
		self.fasta_files_list_to_concat = os.path.join(self.mrbayes_concat_work_dir, "fasta-files-to-concat.txt")
		self.concat_log_filename = os.path.join(self.mrbayes_concat_work_dir,"concat-alignment.out")
		self.concat_log_out_filename = os.path.join(self.mrbayes_concat_work_dir, "concat-alignment.out")
		self.concat_log_err_filename = os.path.join(self.mrbayes_concat_work_dir, "concat-alignment.err")

		self.rewrite_names_dict_filename = os.path.join(self.working_dir,"rewritten-names-concat.json")
		self.names_to_resolved_names = dict()

		self.concat_by_type_workdir = dict()
		self.concat_by_type_workdir[SeqType.Chloroplast] = os.path.join(self.mrbayes_concat_work_dir)
		self.concat_by_type_workdir[SeqType.Mt] = os.path.join(self.mrbayes_concat_work_dir)
		self.concat_by_type_workdir[SeqType.Nuc] = os.path.join(self.mrbayes_concat_work_dir)
		self.concat_by_type_fastas = dict()
		self.concat_by_type_fastas[SeqType.Chloroplast] = os.path.join(self.concat_by_type_workdir[SeqType.Chloroplast] ,"chloroplast-all.fasta")
		self.concat_by_type_fastas[SeqType.Mt] = os.path.join(self.concat_by_type_workdir[SeqType.Mt] ,"Mt-all.fasta")
		self.concat_by_type_fastas[SeqType.Nuc] = os.path.join(self.concat_by_type_workdir[SeqType.Nuc] ,"Nuc-all.fasta")
		self.concat_by_type_outgroup_fasta = dict()
		self.concat_by_type_outgroup_fasta[SeqType.Chloroplast] = os.path.join(self.concat_by_type_workdir[SeqType.Chloroplast] ,"chloroplast-outgroup.fasta")
		self.concat_by_type_outgroup_fasta[SeqType.Mt] = os.path.join(self.concat_by_type_workdir[SeqType.Mt] ,"Mt-outgroup.fasta")
		self.concat_by_type_outgroup_fasta[SeqType.Nuc] = os.path.join(self.concat_by_type_workdir[SeqType.Nuc] ,"Nuc-outgroup.fasta")
		self.concat_by_type_seq_and_outgroup_fasta = dict()
		self.concat_by_type_seq_and_outgroup_fasta[SeqType.Chloroplast] = os.path.join(self.concat_by_type_workdir[SeqType.Chloroplast] ,"chloroplast-seq-and-out.fasta")
		self.concat_by_type_seq_and_outgroup_fasta[SeqType.Mt] = os.path.join(self.concat_by_type_workdir[SeqType.Mt] ,"Mt-seq-and-out.fasta")
		self.concat_by_type_seq_and_outgroup_fasta[SeqType.Nuc] = os.path.join(self.concat_by_type_workdir[SeqType.Nuc] ,"Nuc-seq-and-out.fasta")


		#
		# This part is for the MB part
		#
		self.mb_config_file = os.path.join(self.mrbayes_concat_work_dir, "mb_config.nex")
		self.mb_seq_file = os.path.join(self.mrbayes_concat_work_dir, "mb_final_seq.nex")
		self.mb_out_file = os.path.join(self.mrbayes_concat_work_dir, "mb.out")
		self.mb_out_consensus_file = os.path.join(self.mrbayes_concat_work_dir, "mb.out.con.tre")
		self.mb_out_consensus_file_newick = os.path.join(self.mrbayes_concat_work_dir, "mb.out.con.newick")
		self.mb_out_consensus_file_newick2 = os.path.join(self.mrbayes_concat_work_dir, "mb.out.con.2.newick")

		self.mb_concat_log_out_filename = os.path.join(self.mrbayes_concat_work_dir, "mb-log.out")
		self.mb_concat_log_err_filename = os.path.join(self.mrbayes_concat_work_dir, "mb-log.err")

		self.create_nexus_concat_err_filename = os.path.join(self.mrbayes_concat_work_dir, "mb_concat.err")

		self.mrbayes_concat_rnr_work_dir = os.path.join(self.mrbayes_concat_work_dir,"rnr")
		#
		# This part is for the chromosome numbers part
		#
		self.chr_counts_filename = os.path.join(working_dir,"/chromosome_counts.fa")
		# This part is for the chromevol part
		self.out_trees_filename = os.path.join(working_dir,"parsemb_trees.tre")
		self.out_consensus_tree_filename = os.path.join(working_dir,"parsemb_con_tree.newick")
		self.out_map_tree_filename = os.path.join(working_dir,"parsemb_map_tree.tre")
		self.mr_bayes_post_process_path = os.path.join(working_dir,"mrbayes")

		# ChromEvol
		self.chromevol_workdir = None

		self.pipeline_err_filename = os.path.join(working_dir,"qsub-a_%s.ERR" % id)

	def init_chromevol_dirs(self,chromevol_workdir):
		# ChromEvol
		self.chromevol_workdir = chromevol_workdir
		self.chromevol_job_submit = os.path.join(self.chromevol_workdir,"submit-run-chromevol.sh")
		self.chromevol_control = os.path.join(self.chromevol_workdir,"PIP_control")
		self.chromevol_stopped = os.path.join(self.chromevol_workdir ,"chromevol_stopped")
		self.chromevol_out_dir = os.path.join(self.chromevol_workdir,"chromevol_out")
		self.chromevol_ploidy = os.path.join(self.chromevol_out_dir ,"ploidy.txt")
		self.chromevol_counts = os.path.join(self.chromevol_workdir ,"chromosome-counts.fasta")
		self.chromevol_guided_tree_filename = os.path.join(self.chromevol_workdir ,"guided-tree.fasta")


	def get_parse_filenames(self, base_path):
		return os.path.join(base_path,"parsemb_con_tree.newick"),os.path.join(base_path,"parsemb_trees.tre"), os.path.join(base_path,"parsemb_map_tree.tre"),

	@staticmethod
	def get_mb_filenames(mb_base_dir, mb_basename = "mb.out"):
		mb_files = list()
		mb_files.append(os.path.join(mb_base_dir,"%s.ckp" % mb_basename))
		#mb_files.append(os.path.join(mb_base_dir,"%s.ckp~" % mb_basename))
		mb_files.append(os.path.join(mb_base_dir,"%s.con.tre" % mb_basename))
		mb_files.append(os.path.join(mb_base_dir,"%s.mcmc" % mb_basename))
		mb_files.append(os.path.join(mb_base_dir,"%s.parts" % mb_basename))
		mb_files.append(os.path.join(mb_base_dir,"%s.run1.p" % mb_basename))
		mb_files.append(os.path.join(mb_base_dir,"%s.run1.t"% mb_basename))
		mb_files.append(os.path.join(mb_base_dir,"%s.run2.p" % mb_basename))
		mb_files.append(os.path.join(mb_base_dir,"%s.run2.t" % mb_basename))
		mb_files.append(os.path.join(mb_base_dir,"%s.trprobs" % mb_basename))
		mb_files.append(os.path.join(mb_base_dir,"%s.tstat" % mb_basename))
		mb_files.append(os.path.join(mb_base_dir,"%s.vstat" % mb_basename))
		return mb_files

	@staticmethod
	def get_mr_aic_filenames_in_dir(dir):
		print("Looking for mr aic files in %s" % dir)
		matches = []
		for root, dirnames, filenames in os.walk(dir):
			for filename in fnmatch.filter(filenames, '*phy.phy*'):
				print("in %s" % filename)
				matches.append(os.path.join(root, filename))
		return matches

	def get_mr_aic_filenames_in_all_clusters(self):
		return PathHelper.get_mr_aic_filenames_in_dir(self.working_dir)
