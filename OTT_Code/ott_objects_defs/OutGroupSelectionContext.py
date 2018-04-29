from ott_objects_defs.ploidbCommon import *
#from ott_objects_defs.ploidbCommon import logger, get_cluster_ids_total_data_matrix, get_cluster_ids_total_seq_len

__author__ = 'Moshe'

class OutGroupSelectionContext:
	def __init__(self,ploidb_context,outgroup_name,outgroup_type,all_cids_with_outgroup_seq,all_cids_without_outgroup_seq):
		self.all_cids_with_outgroup_seq = all_cids_with_outgroup_seq.copy()
		self.all_cids_without_outgroup_seq = all_cids_without_outgroup_seq.copy()
		self.outgroup_name = outgroup_name
		self.outgroup_type = outgroup_type
		self.ploidb_context = ploidb_context

		self.candidate_cids_with_outgroup_seq = all_cids_with_outgroup_seq.copy()
		self.candidate_cids_without_outgroup_seq = all_cids_without_outgroup_seq.copy()
		self.not_selected_cids = list()
		self.selected_cids_with_outgroup_seq = list()
		self.selected_cids_without_outgroup_seq = list()
		# This is used in cases where the last cluster that was checked was without outgroup and broke the ratio between the outgroup len and the concat len
		self.next_selection_must_have_outgroup = False


	def __setstate__(self,state):
		self.ploidb_context = None

	def __getstate__(self):
		odict = self.__dict__.copy()
		if 'ploidb_context' in odict:
			del odict['ploidb_context']
		return odict

	def get_selected_cids(self):
		select_cids = self.selected_cids_with_outgroup_seq + self.selected_cids_without_outgroup_seq
		#logger.debug("Returning selected cids %s" % ",".join(select_cids))
		return select_cids

	def add_cluster(self,cid):
		logger.debug("Adding the following cluster for the outgroup list %s" % self.ploidb_context.get_cluster_by_id(cid))
		if cid in self.all_cids_with_outgroup_seq:
			self.selected_cids_with_outgroup_seq.append(cid)
			self.candidate_cids_with_outgroup_seq.remove(cid)
		elif cid in self.all_cids_without_outgroup_seq:
			self.selected_cids_without_outgroup_seq.append(cid)
			self.candidate_cids_without_outgroup_seq.remove(cid)
		else:
			logger.critical("Error cid %s appears in no known list" % cid)
		# Cluster was added -reset this flag
		self.next_selection_must_have_outgroup = False

	def add_cluster_to_not_selected(self,cid):
		self.not_selected_cids.append(cid)
		if cid in self.candidate_cids_without_outgroup_seq:
			self.candidate_cids_without_outgroup_seq.remove(cid)
		if cid in self.candidate_cids_with_outgroup_seq:
			self.candidate_cids_with_outgroup_seq.remove(cid)


	# Determines the score of the outgroup - how "good" is it
	def get_outgroup_score(self):
		# Basic score is based on the data matrix
		score = self.get_selected_cids_data_matrix()

		# Penalty in case the outgroup is genus and not a real organism - it ensure that a "genus" outgroup score is smaller than organism
		if self.outgroup_type == "genus":
			score -= 0.01
		return score

	def get_selected_cids_data_matrix(self):
		return get_cluster_ids_total_data_matrix(self.ploidb_context,self.get_selected_cids())

	def get_outgroup_coverage(self):
		outgroup_coverage_len = get_cluster_ids_total_seq_len(self.ploidb_context,self.selected_cids_with_outgroup_seq)
		total_coverage_len = get_cluster_ids_total_seq_len(self.ploidb_context,self.get_selected_cids())
		if total_coverage_len == 0:
			return 0
		else:
			return (1.0 *outgroup_coverage_len) / total_coverage_len

	def are_candidates_left(self):
		return len(self.candidate_cids_with_outgroup_seq) + len(self.candidate_cids_without_outgroup_seq) > 0

	def get_clusters_with_outgroup(self):
		return self.ploidb_context.get_clusters_by_ids(self.selected_cids_with_outgroup_seq)

	def get_species_in_all_selected_clusters(self):
		species_in_all_selected_clusters = set()

		for cc in self.ploidb_context.get_clusters_by_ids(self.get_selected_cids()):
			species_in_all_selected_clusters.update(cc.get_list_of_species())
		return species_in_all_selected_clusters

	def __str__(self):
		if self.ploidb_context is None:
			return ""
		else:
			with_cid_str = ",".join(self.selected_cids_with_outgroup_seq)
			without_cid_str = ",".join(self.selected_cids_without_outgroup_seq)
			return "[%s (%s)]: score=%d outgroup_coverage=%f data_matrix=%d with_outgroup=%d without=%s with_ids[%s] without_ids[%s]" % \
				   (self.outgroup_name,self.outgroup_type,self.get_outgroup_score(),self.get_outgroup_coverage(),self.get_selected_cids_data_matrix(), len(self.selected_cids_with_outgroup_seq),len(self.selected_cids_without_outgroup_seq),with_cid_str,without_cid_str)
