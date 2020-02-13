__author__ = 'ItayM3'

#--------------------------------------------------------------------------
# (Michal) This version of Name resolution can work with any database
# that has a table with a Name column (that you trust) and thus return
# a matched name, conducting spelling mistakes corrections /
# fuzzy matching and such, using Taxonome
#--------------------------------------------------------------------------

import argparse
import csv
from taxonome.taxa.collection import run_match_taxa
from taxonome import tracker
from taxonome.taxa import file_csv, TaxonSet, Taxon, name_selector
import sqlite3
from taxonome.taxa.file_csv import load_taxa
from taxonome.taxa.name_selector import NameSelector
from taxonome.tracker import CSVTaxaTracker, _flatten_name, CSVListMatches
import pandas as pd

import unicodedata


DEBUG_PRINT = 0

#Dictionaries for Name data:
Glob_Name_TaxId_dict={}
Glob_Name_ParentTaxId_dict={}
Glob_Name_Type_dict={}

__author__ = 'moshe'

class FirstNameSelector(NameSelector):
	def user_select(self, name_options, name, allow_retype, tracker):
		name_to_return = name_options[0]
		try:
			if DEBUG_PRINT == 1: print("handling multiple names for %s - return %s" % (name,name_to_return))
		except:
			if DEBUG_PRINT == 1: print("ofer bug fix 14/02/2016. print caused exception")

		return name_to_return


class MyCSVListMatches(CSVListMatches):
	taxon = None
	matched_name = None

	def __init__(self, fileobj, fieldnames, header=True):
		# Full Name/Authority - accepted name/authority.
		# Name/Authority - accepted name/authority after removing non ascii chars
		# Matched Name - closest sysnoym (or accepted) name found in taxonome DB
		fieldnames = ["TaxId","ParentTaxId","Type","Name", "Authority", "Original name", "Original authority","Coded Name", "Coded Authority", "Score","Matched Name"] + fieldnames
		self.writer = csv.DictWriter(fileobj, fieldnames, extrasaction='ignore')
		if header:
			self.writer.writeheader()

	def start_taxon(self, tax):
		super().start_taxon(tax)
		self.taxon = tax

	def name_transform(self, name, newname, event, **kwargs):
		if DEBUG_PRINT == 1: print(name)
		if DEBUG_PRINT == 1: print(newname)
		if DEBUG_PRINT == 1: print(event)
		super().name_transform(name, newname, event, **kwargs)
		# The idea is: (1) always have a matched_name (2) matched name should be the intermediate step i.e. before synonyms/accepted names
		if (self.matched_name is None) or (event != 'synonymy' and event != 'preferring accepted name'):
			self.matched_name = newname


	def reset(self):
		if self.fromname:
			score = self.fuzzyscore if self.toname else None
			d = dict(self.taxon.info)
			d['Original name'], d['Original authority'] = _flatten_name(self.taxon.name)
			d['Name'], d['Authority']= _flatten_name(self.toname)
			d['Score'] = score


			d['Matched Name'] = None
			matched_name_noauth, matched_author = _flatten_name(self.matched_name)
			#if matched_name_noauth is not None:
			#	d['Matched Name'] = "%s %s" % (matched_name_noauth, matched_author)
			if matched_name_noauth is not None:
				d['Matched Name'] = "%s" % (matched_name_noauth)
				#d['Matched Name'] = "%s" % (self.matched_name)

			for key in Glob_Name_TaxId_dict.keys():
				if DEBUG_PRINT == 1: print(key)
				if DEBUG_PRINT == 1: print(Glob_Name_TaxId_dict[key])#Acanthorrhinum ramosissimum

			if d['Matched Name'] != None and d['Matched Name'] in Glob_Name_TaxId_dict.keys():
				d['TaxId'] = Glob_Name_TaxId_dict[d['Matched Name']]
				d['ParentTaxId'] = Glob_Name_ParentTaxId_dict[d['Matched Name']]
				d['Type'] = Glob_Name_Type_dict[d['Matched Name']]


				d['Coded Name'] = escape_organism_name(d['Name'])
				d['Coded Authority'] = escape_organism_name(d['Authority'])
				# Replace hybrid_marker '×' with '× ' (add the missing space): added by Michal D.
				d['Name'] = handle_Species_hybrid_marker_name(d['Name'])
				d['Original name'] = handle_Species_hybrid_marker_name(d['Original name'])
				d['Matched Name'] = handle_Species_hybrid_marker_name(d['Matched Name'])
			else:
				d['TaxId'] = 'xxx'
				d['ParentTaxId'] = 'xxx'
				d['Type'] = 'xxx'


				d['Coded Name'] = escape_organism_name(d['Name'])
				d['Coded Authority'] = escape_organism_name(d['Authority'])
				# Replace hybrid_marker '×' with '× ' (add the missing space): added by Michal D.
				d['Name'] = 'xxx'
				d['Original name'] = handle_Species_hybrid_marker_name(d['Original name'])
				d['Matched Name'] = None
				d['Score']='0'

			self.writer.writerow(d)

		self.fromname = None
		self.toname = None
		self.taxon = None
		self.matched_name = None


def generate_csv_with_relevant_genus(db_filename, input_filename,out_file,namefield,authfield):
	with open(input_filename, encoding='utf-8', errors='ignore') as f:
		input_taxa = load_taxa(f, namefield=namefield, authfield=authfield)


	#check which DataBase:
	if 'TPL' in db_filename:
		db_type='PLT'
	elif 'COL' in db_filename:
		db_type='COL'
	else:
		db_type='NCBI'

	if DEBUG_PRINT == 1: print(input_filename)
	genus_set = set()
	for name in input_taxa.names:
		if DEBUG_PRINT == 1: print(name)
		species_name_parts = name.split(' ')
		genus_set.add(species_name_parts[0].strip())

	if DEBUG_PRINT == 1: print("Using the following genera: %s " % ",".join(genus_set))
	where_clause = ""

	conn = sqlite3.connect(db_filename)
	curs = conn.cursor()


	with open(out_file,mode="w",encoding='utf8',newline='') as out_handle:
		writer = csv.writer(out_handle,delimiter=',')
		writer.writerow(['Name','TaxId','ParentTaxId','Type','Author','Accepted_name','Accepted_author'])

		add_genus_count = 0
		while len(genus_set) > 0:
			# TODO: put back with lower etc
			genus = genus_set.pop()
			if db_type == 'NCBI':
				where_clause += " OR Name like '%s%%'" % genus
				where_clause += " OR lower(Name) like lower('%%%s%%')" % genus
			elif db_type == 'PLT':  #The plant list
				where_clause += " OR Specie like '%s%%'" % genus
				where_clause += " OR lower(Specie) like lower('%%%s%%')" % genus
			elif db_type == 'COL':  #Catalog of Life
				where_clause += " OR scientificName like '%s%%'" % genus
				where_clause += " OR lower(scientificName) like lower('%%%s%%')" % genus
			add_genus_count += 1

			if add_genus_count > 30 or len(genus_set) == 0:
				if db_type == 'NCBI':
					query_by_genus = "SELECT Name,TaxID,ParentTaxID,Type,'',Name,'' "\
								 "FROM ncbi_names_db "\
								 "WHERE 1=0 %s" % where_clause
				elif db_type == 'PLT':
					query_by_genus = "SELECT Specie,Tpl_ID,Accepted_ID,Status,'',Specie,'' "\
								 "FROM TPL_Name_vs_Accepted "\
								 "WHERE 1=0 %s ORDER BY Status ASC" % where_clause
				elif db_type == 'COL':
					query_by_genus = "SELECT scientificName,taxonID,acceptedNameUsageID,taxonomicStatus,'',scientificName,'' "\
								 "FROM col_names_table "\
								 "WHERE 1=0 %s ORDER BY taxonomicStatus ASC" % where_clause

				curs.execute(query_by_genus)
				rows = curs.fetchall()
				writer.writerows(rows)
				if DEBUG_PRINT == 1: print("%d results were returned running query: %s " % (len(rows),query_by_genus))

				add_genus_count = 0
				where_clause = ""

	conn.close()


def update_dict_glob(synonym_filename):  #Input_NR_check.csv-syn.csv
	#Name,TaxId,ParentTaxId,Type,Author,Accepted_name,Accepted_author
	df = pd.read_csv(synonym_filename, names=['Name','TaxId','ParentTaxId','Type','Author','Accepted_name','Accepted_author'])
	for index, row in df.iterrows():
		name = row['Name'].rstrip() #you can also use df['column_name']
		tpl_id = row['TaxId'] #you can also use df['column_name']
		parent_id = row['ParentTaxId']
		Type_name = row['Type'] #you can also use df['column_name']
		if name not in Glob_Name_TaxId_dict.keys():
			Glob_Name_TaxId_dict[name]=tpl_id
			Glob_Name_ParentTaxId_dict[name]=parent_id
			Glob_Name_Type_dict[name]=Type_name


def write_synonym_results(input_filename, synonym_filename,mappings_file,log_file,namefield='Name', authfield='Author'):
	with open(synonym_filename,encoding='utf8', errors='ignore') as filehandle:
		taxonset = file_csv.load_synonyms(filehandle, synnamefield='Name', synauthfield='Author',
										  accnamefield='Accepted_name', accauthfield='Accepted_author', get_info=True)

	trackers = []
	files = []
	if DEBUG_PRINT == 1: print("Openning %s" % mappings_file)
	f = open(mappings_file, "w", encoding='utf-8', newline='')
	files.append(f)
	trackers.append(MyCSVListMatches(f,["Id"]))

	if DEBUG_PRINT == 1: print("Openning %s" % log_file)
	f = open(log_file, "w", encoding='utf-8', newline='')
	files.append(f)
	trackers.append(tracker.CSVTracker(f))

	if DEBUG_PRINT == 1: print("Loading taxa from input file %s" % input_filename)
	with open(input_filename, encoding='utf-8', errors='ignore') as f:
		input_taxa = load_taxa(f, namefield=namefield, authfield=authfield)

	if DEBUG_PRINT == 1: print("Running match taxa. inpout taxa len=%d synonym len=%d" % (len(input_taxa),len(taxonset)))
	#run_match_taxa (input_taxa,taxonset, tracker=trackers,nameselector=FirstNameSelector(),prefer_accepted='all')
	run_match_taxa (input_taxa,taxonset, tracker=trackers,nameselector=FirstNameSelector(),prefer_accepted='noauth')

	for f in files: f.close()


def handle_Species_hybrid_marker_name(organism_name):
	if organism_name is None:
		return None

	escaped_organism = organism_name
	escaped_organism = escaped_organism.replace("×", "× ")

	return escaped_organism


def escape_organism_name(organism_name):
	if organism_name is None:
		return None


	escaped_organism = organism_name

	escaped_organism = escaped_organism.replace("×", "x ")
	escaped_organism = escaped_organism.replace(" ", "_")
	escaped_organism = escaped_organism.replace("]", "")
	escaped_organism = escaped_organism.replace("[", "")
	escaped_organism = escaped_organism.replace(",", "_")
	escaped_organism = escaped_organism.replace("-", "_")
	escaped_organism = escaped_organism.replace("'", "_")
	escaped_organism = escaped_organism.replace("/", "_")
	escaped_organism = escaped_organism.replace("&", "AND")
	escaped_organism = unicodedata.normalize('NFKD', escaped_organism).encode('ascii','ignore')
	escaped_organism = str(escaped_organism, encoding='ascii')

	return escaped_organism

def do_resolve_names(db_filename,input_filename,results_filename,log_filename,namefield,authfield):
	genera_synonyms_csv_filename = input_filename + "-syn.csv"
	generate_csv_with_relevant_genus(db_filename=db_filename, input_filename=input_filename,out_file=genera_synonyms_csv_filename,namefield=namefield,authfield=authfield)
	update_dict_glob(genera_synonyms_csv_filename)
	write_synonym_results(input_filename=input_filename,synonym_filename = genera_synonyms_csv_filename ,mappings_file = results_filename,log_file = log_filename,namefield=namefield,authfield=authfield)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Main ploiDB script - builds phylogenetic trees per genus')
	parser.add_argument('--db-filename','-d', help='SQLite DB filename containing the names from plant list', required=True)
	parser.add_argument('--input-filename','-i', help='File with the list of genera to run the matches on', required=True, default=None)
	parser.add_argument('--results-filename','-r', help='Output filename for the results', required=True, default=None)
	parser.add_argument('--log-filename','-l', help='log filename', required=True, default=None)
	parser.add_argument('--authfield','-a', help='author field. Can be name of field, True or None. \
				True means that the author is part of the name, None means there is no author and other string indicates the name of the column in the csv with the author name', required=False, default=None)

	args = parser.parse_args()

	authfield = args.authfield
	if authfield == 'None' or authfield is None:
		if DEBUG_PRINT == 1: print ("Author field is %s. Author name is not provided at all" % authfield)
		authfield = None
	elif authfield == 'True' or authfield is True:
		if DEBUG_PRINT == 1: print ("Author field is %s. Author name is part of the species name" % authfield)
		authfield = True
	elif authfield == 'False' or authfield is False:
		if DEBUG_PRINT == 1: print ("Author field is %s. Not sure what it means... Don't use this False" % authfield)
		authfield = False
	else:
		if DEBUG_PRINT == 1: print ("Author field is %s - author name will be in this column in the csv" % authfield)

	do_resolve_names(args.db_filename,args.input_filename,args.results_filename,args.log_filename,'Name',authfield)

