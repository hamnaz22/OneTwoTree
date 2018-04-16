from __future__ import division
import csv
import logging
import sys
import math

__author__ = 'Shiran'


def get_avg(l):
	avg = sum(l)/len(l)
	return avg


def get_var(l):
	avg = get_avg(l)
	dist_from_avg = list(map(lambda x: (x - avg)**2, l))
	var = sum(dist_from_avg)/len(dist_from_avg)
	return var


def get_std(l):
	var = get_var(l)
	std = math.sqrt(var)
	return std


def sum_of_squares(l):
	return sum([c**2 for c in l])


def calc_division(numerator, denominator):
	TOL = 0.0000000001
	if abs(denominator) > TOL and abs(numerator) > TOL:
		return numerator/denominator
	elif (abs(denominator) < TOL) and (abs(numerator) > TOL):
		sign = numerator/abs(numerator)
		return sign*10**6
	else:
		return 0


def init_commandline_logger(logger):
	logger.setLevel(logging.DEBUG)
	# create console handler and set level to debug
	ch = logging.StreamHandler(sys.stdout)
	ch.setLevel(logging.DEBUG)
	# create formatter
	formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	# add formatter to ch
	ch.setFormatter(formatter)
	# add ch to logger
	logger.addHandler(ch)


def export_dict_to_csv(dict, csv_name, headings):
	with open(csv_name, "w", newline='') as fp:
		w = csv.writer(fp, dialect='excel', quoting=csv.QUOTE_NONNUMERIC)
		w.writerow(headings)
		for key in dict.keys():
			w.writerow([key] + dict[key])


def import_csv_to_dict(csv_name):
	dict ={}
	with open(csv_name, "r") as fp:
		r = csv.reader(fp)
		for x in r:
			dict[x[0]] = x[1:]
	dict.pop("cluster_ID")
	return dict