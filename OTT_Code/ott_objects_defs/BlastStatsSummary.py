__author__ = 'moshe'

"""
Class for holding summary for blast execution in order to find an outgroup
"""

class BlastStatsSummary:

	def __init__(self):
		self.desc = None
		self.number_of_results = 0

		self.min_eval = -1.0
		self.max_eval = -1.0
		self.avg_eval = -1.0
		self.percentile_eval = list()
		self.percentile_eval_str = ""
		self.low_percentile_eval = -1.0
		self.high_percentile_eval = -1.0
		self.stddev_eval = -1.0

		self.min_score = -1.0
		self.max_score = -1.0
		self.avg_score = -1.0
		self.percentile_score = list()
		self.percentile_score_str = ""
		self.low_percentile_score = -1.0
		self.high_percentile_score = -1.0
		self.stddev_score = None

	def __str__(self):
		return "Stats for [%s] [%d results] " \
			   "[eval: min=%f max=%f avg=%f low_percentile_score=%f high_percentile_score=%f stddev=%f] " \
			   "[score: min=%f max=%f avg=%f low_percentile_score=%f high_percentile_score=%f stddev=%f]" % (self.desc, self.number_of_results, self.min_eval,
						self.max_eval, self.avg_eval,self.low_percentile_eval,self.high_percentile_eval,self.stddev_eval,
						self.min_score, self.max_score, self.avg_score, self.low_percentile_score,self.high_percentile_score,self.stddev_score)

