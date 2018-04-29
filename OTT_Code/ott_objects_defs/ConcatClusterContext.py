from ott_objects_defs.ClusterContext import * #ClusterContext
from ott_objects_defs.ploidbCommon import * #concat_fasta_cluster
from Bio import SeqIO

__author__ = 'moshe'

"""
This class is a data structure for holding a list of clusters. Currently it is used for holding all the chloroplast clusters
in a single data structure since we want to treat them eventually (in MrBayes) as a single cluster
"""


class ConcatClusterContext(ClusterContext):
	CHLOROPLAST_CLUSTER_INDEX =  999
	MITOCHONDRIAL_CLUSTER_INDEX =  666
	NUCLEUS_CLUSTER_INDEX =  333

	def __init__(self,cluster_index,all_seqs_fasta_filename,genus_cluster,clusters_list,cluster_type):
		super(ConcatClusterContext,self).__init__(cluster_index,all_seqs_fasta_filename,genus_cluster)
		self.cluster_list = clusters_list
		self.cluster_type = cluster_type
		self.gene_name = None

	def init_cluster_id(self):
		if self.cluster_type == SeqType.Chloroplast:
			self.cluster_id = self.CHLOROPLAST_CLUSTER_INDEX
		elif self.cluster_type == SeqType.Mt:
			self.cluster_id = self.MITOCHONDRIAL_CLUSTER_INDEX
		if self.cluster_type == SeqType.Nuc:
			self.cluster_id = self.NUCLEUS_CLUSTER_INDEX

	def init_cluster_type(self):
		# was determined in the creation
		pass

	def init_cluster(self,genus_context):
		concat_fasta_cluster(clusters=self.cluster_list,workdir=self.cluter_work_dir,concat_fasta=self.original_fasta_filename,
							 concat_outgroup_fasta=self.outgroup_seq_filename)
		super(ConcatClusterContext, self).init_cluster(genus_context)
		# Setting this property according to the first one in the list so that later on, when its used, there won't be a problem
		self.no_of_outgroup_sequence_for_concat =self.cluster_list[0].no_of_outgroup_sequence_for_concat
		self.outgroup_sequence_for_concat_short_desc_filename = self.outgroup_seq_filename
		self.all_seqs_for_concat_tree_creation_fasta_filename = self.original_fasta_filename
		self.init_cluster_id()

	def init_blast(self,genus_context):
		# Nothing to do - blast was done for each cluster separately
		pass

	def get_number_of_different_species(self):
		if self.number_of_different_species is None:
			all_seqs = list(SeqIO.parse(self.original_fasta_filename, "fasta"))
			self.number_of_different_species = len(all_seqs)
		return self.number_of_different_species


	def get_data_matrix_size(self,estimated=False):
		total_data_matrix_size = 0
		for cc in self.cluster_list:
			total_data_matrix_size += cc.get_data_matrix_size(estimated)
		return total_data_matrix_size

	def get_cluster_length(self,estimated=False):
		total_cluster_length = 0
		for cc in self.cluster_list:
			total_cluster_length += cc.get_cluster_length(estimated)
		return total_cluster_length


