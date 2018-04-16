# handle ITS sequences for a given genus
use strict;
use Bio::Index::GenBank;
use FindBin qw($Bin);
use lib $Bin; # will have the same path where the .pl is
use ITS;
use Guidance;

if(@ARGV < 7){
	die "usage:
	(1) output directory of a single genus
	(2) sequences in fasta format
	(3) genBank file
	(4) scripts dir - location of the scripts used by this one
	(5) config .ini file
	(6) guidance flag
	(7) msa_software
	";
}

my ($genusOutDir, $fasta, $GB, $scriptsDir, $configFile, $guidanceFlag, $msa_software) = @ARGV;

my $gbIndex = Bio::Index::GenBank->new (-filename => $GB . ".idx", -write_flag => 1);
$gbIndex->make_index($GB);

mkdir $genusOutDir unless (-d $genusOutDir);

#my ($log, $genus) = ($genusOutDir =~ m/^(.+)\/([^\/]+)$/);
my $pickFromFastalog = "$genusOutDir/pickOneSeqFromFasta.log";
my $pickFromMSAlog = "$genusOutDir/pickOneSeqFromMSA.log";
my $its1fasta = "$genusOutDir/ITS1_only.fasta"; my $its2fasta = "$genusOutDir/ITS2_only.fasta"; my $itscombfasta = "$genusOutDir/combined.fasta";
# Filtering the sequences so that for each species we would have at most one seq of each ITS type (ITS1/ITS2/combined)

print "Filtering $fasta - each species will have at most one ITS seq of each type (ITS1/ITS2/combined) \n";
ITS::pickOneITSTypePerSpeciesFasta($fasta, $gbIndex, $genusOutDir, "$genusOutDir/oneITSTypePerSpecies.fasta", $scriptsDir, $pickFromFastalog);

# separate the sequences, after they were filtered, to separate files that contain ITS1, ITS2 and ITS1+ITS2:
print "Splitting filtered results by ITS type\n";
my ($ITS1count, $ITS2count, $combinedCount) = ITS::splitITS("$genusOutDir/oneITSTypePerSpecies.fasta",$its1fasta , $its2fasta,$itscombfasta , $gbIndex);

# if there are no combined sequences in the genus, and the separated ITS1 and ITS2 sequences aren't aligned to each other,
# so sequences from the same species should be appended and NOT merged
my $append = 1;
if ($combinedCount > 0) {
    $append = 0; # merge
}

if($ITS1count > 0 || $ITS2count > 0 || $combinedCount > 0){
    # Do MSA for the ITS seqs
    my $MSAfileName;
    print "Found ITS seqs - starting MSA for ITS\n";
	#'ClustalOmega' or 'MAFFT'
	if ($msa_software eq 'ClustalOmega'){
	    print "MSA software -> ClustalOmega\n";
		$MSAfileName = ITS::ITS_CLUSTALO($genusOutDir, $ITS1count, $ITS2count, $combinedCount, $its1fasta, $its2fasta, $itscombfasta,$scriptsDir);
	}else{
	    print "MSA software -> MAFFT\n";
		$MSAfileName = ITS::ITS_MSA($genusOutDir, $ITS1count, $ITS2count, $combinedCount, $its1fasta, $its2fasta, $itscombfasta,$scriptsDir);
	}
	print "After MSA - selecting a single ITS seq per species. MSA file: $MSAfileName\n";
    ITS::pickOneSeqPerSpeciesMSA($MSAfileName, $gbIndex, $genusOutDir,  "$genusOutDir/oneSeqPerSpecies.msa", $scriptsDir, $pickFromMSAlog,$append);

	# Michal: NOT sure why we need this code? !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if($MSAfileName =~ m/combined\+sep\.msa/){

		print "About to count number of seqs in $genusOutDir/SEP_ITS1+ITS2.fasta\n";
		my $count_sep=0;
		my $in=Bio::SeqIO->new(-format=>'fasta', -file=>"$genusOutDir/SEP_ITS1+ITS2.fasta");
		while(my $seq=$in->next_seq){
			$count_sep++;
		}

		print "Found $count_sep seqs in $genusOutDir/SEP_ITS1+ITS2.fasta\n";
		my $count_comb=0;
		my $in=Bio::SeqIO->new(-format=>'fasta', -file=>$itscombfasta);
		while(my $seq=$in->next_seq){
			$count_comb++;
		}
		print "Found $count_comb seqs in $itscombfasta\n";
		print "GUIDANCE Flag is set to: $guidanceFlag\n";
		if ($guidanceFlag eq "True"){
			if ($count_comb < 5 || $count_sep < 5) {
				print "NOT Calling GUIDANCE since there are only $count_sep and $count_comb sequences\n";
			} else {


				# TODO: calling GUIDANCE should be done with the correct files...

				Guidance::runGuidance($itscombfasta, "$genusOutDir/SEP_ITS1+ITS2.fasta", "$genusOutDir/GuidnaceOutput", $configFile);
				Guidance::checkGuidance($MSAfileName, $MSAfileName.".filtered_out", "$genusOutDir/GuidnaceOutput");

				# TODO - use GUIDANCE results by overriding the MSA file with the filtered results
				#$MSAfileName = $MSAfileName.".filtered_out";
			}
		} else {
			print " Guidance Flag is set to False !!!"
		}

	} else {
	    # TODO - execute GUIDACE even in this case - it should be standatd GUIDANCE execution and not special like the above
	    print "NOT Calling GUIDANCE since its not combination of combined + sep\n";
	}
}