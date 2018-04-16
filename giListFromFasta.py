from ploidbCommon import  *

__author__ = 'moshe'

#
# given a fasta file - write the GI to an output file
#
def giListFromFastaFile(fastaFilename, giListOutFilename):
	logger.info("about to get GI list from " + fastaFilename)
	giListOutHandle = open(giListOutFilename, 'w')
	fastaHandle = open(fastaFilename, "r")
	for seq_record in SeqIO.parse(fastaHandle, "fasta"):
		gi = getPropertyFromFastaSeqHeader(seq_record.description, "gi")
		giListOutHandle.write(gi + "\n")

	giListOutHandle.close()

#
# given a fasta file - return all the GIs in the fasta separated by comma
# michal: added ' around the gi since it has a structure of 'ZX123123123.2' (after replacing the gi number with the accession id)
def gi_list_str_from_fasta_file(fastaFilename):
	logger.info("about to get GI list from " + fastaFilename)
	fastaHandle = open(fastaFilename, "r")
	gis_str = ""
	for seq_record in SeqIO.parse(fastaHandle, "fasta"):
		gi = getPropertyFromFastaSeqHeader(seq_record.description, "gi")
		gis_str += ", " + "'" + gi + "'"

	return gis_str[2:]



