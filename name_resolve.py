import csv
from taxonome.taxa.collection import run_match_taxa
from taxonome import tracker
from taxonome.taxa import file_csv, TaxonSet, Taxon, name_selector
import sqlite3
from taxonome.taxa.file_csv import load_taxa
from taxonome.taxa.name_selector import NameSelector
from taxonome.tracker import CSVTaxaTracker, _flatten_name, CSVListMatches
import unicodedata


"""
Given a csv - execute name resultion using Taxonome
"""

__author__ = 'moshe'


class FirstNameSelector(NameSelector):
	def user_select(self, name_options, name, allow_retype, tracker):
		name_to_return = name_options[0]
		#print("handling multiple names for %s - return %s" % (name,name_to_return))
		return name_to_return


class MyCSVListMatches(CSVListMatches):
	taxon = None
	matched_name = None

	def __init__(self, fileobj, fieldnames, header=True):
		# Full Name/Authority - accepted name/authority.
		# Name/Authority - accepted name/authority after removing non ascii chars
		# Matched Name - closest sysnoym (or accepted) name found in taxonome DB
		fieldnames = ["Name", "Authority", "Original name", "Original authority","Coded Name", "Coded Authority", "Score","Matched Name"] + fieldnames
		self.writer = csv.DictWriter(fileobj, fieldnames, extrasaction='ignore')
		if header:
			self.writer.writeheader()

	def start_taxon(self, tax):
		super().start_taxon(tax)
		self.taxon = tax

	def name_transform(self, name, newname, event, **kwargs):
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
			if matched_name_noauth is not None:
				d['Matched Name'] = "%s %s" % (matched_name_noauth, matched_author)

			d['Coded Name'] = escape_organism_name(d['Name'])
			d['Coded Authority'] = escape_organism_name(d['Authority'])
			self.writer.writerow(d)

		self.fromname = None
		self.toname = None
		self.taxon = None
		self.matched_name = None

"""
	Fetch tbe  names from the DB and write them to a CSV so that taxonome would be able to use it (taxonome only works with csv currently)
"""
def generate_csv_with_relevant_genus(db_filename, input_filename,out_file,namefield,authfield):
	with open(input_filename, encoding='utf-8', errors='ignore') as f:
		input_taxa = load_taxa(f, namefield=namefield, authfield=authfield)

	genus_set = set()
	for name in input_taxa.names:
		species_name_parts = name.split(' ')
		genus_set.add(species_name_parts[0].strip())

	print("Using the following genera: %s " % ",".join(genus_set))
	where_clause = ""

	conn = sqlite3.connect(db_filename)

	curs = conn.cursor()

	with open(out_file,mode="w",encoding='utf8',newline='') as out_handle:
		writer = csv.writer(out_handle,delimiter=',')
		writer.writerow(['Name','Author','Accepted_name','Accepted_author'])

		add_genus_count = 0
		while len(genus_set) > 0:
			# TODO: put back with lower etc
			genus = genus_set.pop()
			where_clause += " OR Name like '%s%%'" % genus
			#where_clause += " OR lower(Name) like lower('%%%s%%')" % genus
			add_genus_count += 1

			if add_genus_count > 30 or len(genus_set) == 0:
				query_by_genus = "SELECT Name, Author,Accepted_name, Accepted_author "\
								 "FROM TPL_names "\
								 "WHERE 1=0 %s" % where_clause

				curs.execute(query_by_genus)
				rows = curs.fetchall()
				writer.writerows(rows)
				print("%d results were returned running query: %s " % (len(rows),query_by_genus))

				add_genus_count = 0
				where_clause = ""

	conn.close()


def write_synonym_results(input_filename, synonym_filename,mappings_file,log_file,namefield='Name', authfield='Author'):
	with open(synonym_filename,encoding='utf8', errors='ignore') as filehandle:
		taxonset = file_csv.load_synonyms(filehandle, synnamefield='Name', synauthfield='Author',
										  accnamefield='Accepted_name', accauthfield='Accepted_author', get_info=True)

	trackers = []
	files = []
	print("Openning %s" % mappings_file)
	f = open(mappings_file, "w", encoding='utf-8', newline='')
	files.append(f)
	trackers.append(MyCSVListMatches(f,["Id"]))

	print("Openning %s" % log_file)
	f = open(log_file, "w", encoding='utf-8', newline='')
	files.append(f)
	trackers.append(tracker.CSVTracker(f))

	print("Loading taxa from input file %s" % input_filename)
	with open(input_filename, encoding='utf-8', errors='ignore') as f:
		input_taxa = load_taxa(f, namefield=namefield, authfield=authfield)

	print("Running match taxa. inpout taxa len=%d synonym len=%d" % (len(input_taxa),len(taxonset)))
	run_match_taxa (input_taxa,taxonset, tracker=trackers,nameselector=FirstNameSelector(),prefer_accepted='all')

	for f in files: f.close()


def escape_organism_name(organism_name):
	if organism_name is None:
		return None

	escaped_organism = organism_name

	#escaped_organism = escaped_organism.replace("×", "X")
	escaped_organism = escaped_organism.replace("×", "x ")
	escaped_organism = escaped_organism.replace(" ", "_")
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
	write_synonym_results(input_filename=input_filename,synonym_filename = genera_synonyms_csv_filename ,mappings_file = results_filename,log_file = log_filename,namefield=namefield,authfield=authfield)
