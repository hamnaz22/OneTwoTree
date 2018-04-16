from BlastStatsSummary import BlastStatsSummary

__author__ = 'moshe'


class OutGroupSelection:

	def __init__(self):
		# num er of blast results that are outside the cluster and not part of the genus
		self.blast_results_outside_cluster = 0
		self.outgroup_found = False

		self.outgroup_organism = None
		self.outgroup_selection_description = None
		self.eval = -1.0
		self.minscore = -1.0

		self.cluster_stats = BlastStatsSummary()
		self.genus_stats = BlastStatsSummary()
		self.non_genus_stats = BlastStatsSummary()