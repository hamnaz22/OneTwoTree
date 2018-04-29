import collections
import inspect
from json import JSONEncoder
from ott_objects_defs.BlastStatsSummary import BlastStatsSummary
from ott_objects_defs.OutGroupSelectionContext import OutGroupSelectionContext
from ott_objects_defs.ClusterContext import *
from ott_objects_defs.PloiDbContext import *
from ott_objects_defs.ploidbCommon import logger

from ott_objects_defs.OutGroupSelection import OutGroupSelection


__author__ = 'moshe'

class ObjectJSONEncoder(JSONEncoder):
	PloiDbContextAttrs = list()
	PloiDbContextAttrs.append('cluster_contexts')
	PloiDbContextAttrs.append('its_final_cluster')
	PloiDbContextAttrs.append('cluster_contexts_for_loci_trees')
	PloiDbContextAttrs.append('cluster_contexts_for_concat_tree')
	PloiDbContextAttrs.append('number_of_clusters')
	PloiDbContextAttrs.append('mb_out_consensus_file')
	PloiDbContextAttrs.append('concat_seqs_fasta_filename')
	PloiDbContextAttrs.append('chloroplast_cluster')
	PloiDbContextAttrs.append('mitochondrial_cluster')
	PloiDbContextAttrs.append('nucleus_cluster')
	PloiDbContextAttrs.append('names_to_resolved_names')
	PloiDbContextAttrs.append('unresolved_species_names')
	PloiDbContextAttrs.append('rnr_concat_filtered_species')
	PloiDbContextAttrs.append('cluster_contexts_for_concat_tree_candidates')
	PloiDbContextAttrs.append('top_outgroup_contexts')
	PloiDbContextAttrs.append('selected_outgroup')
	PloiDbContextAttrs.append('full_cluster_contexts')
	PloiDbContextAttrs.append('irrelevant_cluster_contexts')

	ClusterContextBaseAttrs = list()
	ClusterContextBaseAttrs.append('genus_working_dir')
	ClusterContextBaseAttrs.append('genus_name')
	ClusterContextBaseAttrs.append('cluster_id')
	ClusterContextBaseAttrs.append('cluster_desc')
	ClusterContextBaseAttrs.append('number_of_different_species')
	ClusterContextBaseAttrs.append('is_its_cluster')
	ClusterContextBaseAttrs.append('cluster_desc')
	ClusterContextBaseAttrs.append('cluster_type')
	ClusterContextBaseAttrs.append('gene_name')
	ClusterContextBaseAttrs.append('seq_length')
	ClusterContextBaseAttrs.append('estimated_seq_length')
	ClusterContextBaseAttrs.append('data_matrix_size')
	ClusterContextBaseAttrs.append('index')
	ClusterContextBaseAttrs.append('no_of_outgroup_sequence_for_concat')
	ClusterContextBaseAttrs.append('outgroup_sequence_for_concat_short_desc_filename')

	ClusterContextFullAttrs = list()
	ClusterContextFullAttrs.append('all_seqs_with_outgroup_filename')
	ClusterContextFullAttrs.append('locus_tree_outgroup_selection')
	ClusterContextFullAttrs.append('all_seqs_for_tree_creation_fasta_filename')
	ClusterContextFullAttrs.append('all_seqs_for_concat_tree_creation_fasta_filename')
	ClusterContextFullAttrs.append('mb_out_consensus_file')
	ClusterContextFullAttrs.append('mb_out_consensus_file')

	def default(self, reject):
		is_not_method = lambda o: not inspect.isroutine(o)
		non_methods = inspect.getmembers(reject, is_not_method)
		non_methods_dict = dict()
		for attr, value in non_methods:
			if not attr.startswith('__'):
				if value is not None and (attr == 'cluster_contexts_for_concat_tree_candidates' or attr.startswith('cluster_contexts_for_concat_tree')
							or attr == 'cluster_contexts_for_loci_trees' or attr == 'full_cluster_contexts' or attr == 'irrelevant_cluster_contexts'):
					logger.debug("FOUND %s" % attr)
					non_methods_dict[attr] = list(str(cc) for cc in value)

				elif attr == 'top_outgroup_contexts':
					non_methods_dict[attr] = list(str(outgroup_context) for outgroup_context in value)

				#elif isinstance(reject, PloiDbContext.PloiDbContext):
				elif isinstance(reject, PloiDbContext):
					if attr in self.PloiDbContextAttrs:
						non_methods_dict[attr] = value
				#elif isinstance(reject, ClusterContext.ClusterContext):
				elif isinstance(reject, ClusterContext):
					if attr in self.ClusterContextBaseAttrs or attr in self.ClusterContextFullAttrs:
						non_methods_dict[attr] = value
				elif isinstance(reject, OutGroupSelection):
					non_methods_dict[attr] = value
				elif isinstance(reject, OutGroupSelectionContext):
					non_methods_dict[attr] = str(value)
				elif isinstance(reject, BlastStatsSummary):
					if not type(value) is list:
						non_methods_dict[attr] = value
				else:
					logger.error("attr %s is not in any known object. reject is %s" % (attr,str(reject)))

		return non_methods_dict
