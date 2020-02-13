
use strict;
use Config::Tiny;

#don't forget: clean the database before
#user guide: http://orthomcl.org/common/downloads/software/v2.0/UserGuide.txt

if ( @ARGV < 9 ) {
	die "usage:
	(1) fasta file contains all the sequences,
	(2) tables names suffix (unique),
	(3) percent match cutoff (for blast) - default: 50,
	(4) inflation parameter (for mcl) - default: 1.5
	(5) directory contains all the scripts
	(6) GenBank file contains the sequences of this genus
	(7) output directory
	(8) config .ini file
	(9) initial filter of blast results - default is 0.5
	(10) orthomcl bin directory
	(11) input taxa num - for filtering clusters
";
}

my ($fastaAllSequences, $tablesNamesSuffix, $percentMatchCutoff, $inflationParam, $scriptsDir, $gbFile, $outDir, $configFile, $filterBlastRatio, $orthomclBinDir, $inputTaxaNum) = @ARGV;    #default: $percentMatchCutoff = 50, $inflationParam = 1.5
my (@seqsByOrganismFiles, $taxonCode, $numOfSequnces, $genusName);

# Create a config
my $Config = Config::Tiny->new;
# Open the config
$Config = Config::Tiny->read( $configFile );
# Reading properties
my $rootproperty = $Config->{_}->{rootproperty};
my $mysqlHostName = $Config->{ott_mysql}->{hostname};
my $mysqlUserName = $Config->{ott_mysql}->{username};
my $mysqlPassword = $Config->{ott_mysql}->{password};

my $mysqlLibPath = $Config->{orthomcl}->{ortho_perl_path};


if ($fastaAllSequences =~ m/([^\/]+)-allseq.fasta/){
	$genusName = $1;
}
else{
	$genusName = "NN";
}

#step 4: orthomclInstallSchema
step4();

#step 5: orthomclAdjustFasta
step5();

#step 6: orthomclFilterFasta
step6();

#Step 7: All-v-all BLAST
step7();

filterBlastOutput();

#step 8: orthomclBlastParser
step8();

#step 9: orthomclLoadBlast
step9();

#step 10: orthomclPairs
#my $ret_val = "yes";
#$ret_val = &step10();
step10();

#step11: orthomclDumpPairsFiles
step11();

#NEW - added by michal for cases that orthomcl failes to find pairs:
my $empty_input = "yes";
my $file_mcl_input = "$outDir/mclInput";
#check if mclInput file is empty and create empty groups file:
if(-z $file_mcl_input){
	print("File mclInput at $outDir is empty, create empty groups file");
	my $groups_empty_file = "$outDir/groups.txt";
    open my $fh,'>', $groups_empty_file or die "Can't open the groups file: $!";
	close($fh);
} else {
	$empty_input = "no";
}

if( $empty_input eq "yes"){

	# parse the output:

	my $onlyDesc = $outDir."_onlyDesc";
	my $seqs = $outDir."_seqs";

	system "perl -w $scriptsDir/getSequencesFromClusters.pl $outDir/groups.txt $seqs $fastaAllSequences $inputTaxaNum"; #get the sequences

	system "perl -w $scriptsDir/getOnlyDesc.pl $seqs $onlyDesc"; #get descriptions of each cluster

	system "perl -w $scriptsDir/collapseCluster.pl $onlyDesc $gbFile $seqs $fastaAllSequences n"; #collapse the clusters

	system "perl -w $scriptsDir/appendClusters.pl $onlyDesc"; #append the collapsed clusters to one file - easier to read
	return

} else {


	#Step 12: mcl
	step12();

	#Step 13: orthomclMclToGroups
	step13();


	# parse the output:

	my $onlyDesc = $outDir."_onlyDesc";
	my $seqs = $outDir."_seqs";

	system "perl -w $scriptsDir/getSequencesFromClusters.pl $outDir/groups.txt $seqs $fastaAllSequences $inputTaxaNum"; #get the sequences

	system "perl -w $scriptsDir/getOnlyDesc.pl $seqs $onlyDesc"; #get descriptions of each cluster

	system "perl -w $scriptsDir/collapseCluster.pl $onlyDesc $gbFile $seqs $fastaAllSequences n"; #collapse the clusters

	system "perl -w $scriptsDir/appendClusters.pl $onlyDesc"; #append the collapsed clusters to one file - easier to read
};
 
 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#						Sub-rotines
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sub step4 {
    print "In step 4\n";

	open( CONFIG4, ">$outDir/config_step4" ) or die "can't open file $outDir/config_step4";

	print CONFIG4 "dbVendor=mysql
dbConnectString=dbi:mysql:mysql_local_infile=1:o_" . $tablesNamesSuffix . ":" . $mysqlHostName . "
dbLogin=". $mysqlUserName. "
dbPassword=" . $mysqlPassword ."
similarSequencesTable=SimilarSequences" . $tablesNamesSuffix . "
orthologTable=Ortholog" . $tablesNamesSuffix . "
inParalogTable=InParalog" . $tablesNamesSuffix . "
coOrthologTable=CoOrtholog" . $tablesNamesSuffix . "
interTaxonMatchView=InterTaxonMatch" . $tablesNamesSuffix . "
percentMatchCutoff=" . $percentMatchCutoff . "
evalueExponentCutoff=-5
oracleIndexTblSpc=NONE";

	close(CONFIG4);
    print "perl -w $orthomclBinDir/orthomclInstallSchema $outDir/config_step4 $outDir/log\n";
	system "perl -w $orthomclBinDir/orthomclInstallSchema $outDir/config_step4 $outDir/log";
}

#==============================================
sub step5 {
    print "In step 5\n";

	mkdir("$outDir/seqsByOrganism");
	system "perl -w $scriptsDir/SplitFastaSeqsBySpecies.pl $fastaAllSequences $outDir/seqsByOrganism"; #split the input fasta file by species

	opendir( SEQS_DIR, "$outDir/seqsByOrganism" ) or die "can't open directory $outDir/seqsByOrganism";

	mkdir("$outDir/compliantFasta");

	@seqsByOrganismFiles = grep { $_ ne '.' && $_ ne '..' } readdir(SEQS_DIR);

	foreach my $organismFile (@seqsByOrganismFiles) {
		if ($organismFile =~ m/(\d+)\.fasta/){
			$taxonCode = $1;
		}
#		print "$taxonCode\n";
		system "perl -w $scriptsDir/orthomclAdjustFasta.pl $taxonCode $outDir/seqsByOrganism/$organismFile 2 $outDir/compliantFasta";
	}
}

#==============================================
sub step6 {
    print "In step 6\n";

	system "perl -w $scriptsDir/orthomclFilterFasta.pl $outDir/compliantFasta 10 20 $outDir";
	#arguments: 10=minimum allowed length of proteins, 20=maximum percent stop codons
}

#==============================================
sub step7 {
    print "In step 7\n";

	open( ALL_SEQS, "<$fastaAllSequences" )
	  or die "can't open file $fastaAllSequences";
	$numOfSequnces = 0;
	while (<ALL_SEQS>) {
		if ( substr( $_, 0, 1 ) eq ">" ) {
			$numOfSequnces++;
		}
	}
	close(ALL_SEQS);

	system "formatdb -i $outDir/goodProteins.fasta -pF -l $outDir/formatdb.log";    #-pF -- nucleotides
	system "blastall -p blastn -d $outDir/goodProteins.fasta -i $outDir/goodProteins.fasta -v 100000 -b 100000 -z $numOfSequnces -e 1e-5 -m 8 > $outDir/blast_all-v-all_total.blastn";

}

#==============================================
sub filterBlastOutput{
    print "Filtering short BLAST results\n";
    # Preprocessing the blast results in order to filter sequences that are not covered by a certain threshold in the blast results
    #system "python $scriptsDir/filterBlastResultsByLength.py -i $outDir/blast_all-v-all_total.blastn -f $outDir/goodProteins.fasta -o $outDir/blast_all-v-all_output.blastn -fbr $filterBlastRatio > $outDir/filterblast.out 2>&1" ;
    system "python $scriptsDir/filterBlastResultsByLength.py -i $outDir/blast_all-v-all_total.blastn -f $outDir/goodProteins.fasta -o $outDir/blast_all-v-all_output.blastn -fbr $filterBlastRatio" ;
	print("python $scriptsDir/filterBlastResultsByLength.py -i $outDir/blast_all-v-all_total.blastn -f $outDir/goodProteins.fasta -o $outDir/blast_all-v-all_output.blastn  -fbr $filterBlastRatio> $outDir/filterblast.out 2>&1\n");
}


sub step8{
    print "In step 8\n";

	my $orthomclBlastParser = $Config->{orthomcl}->{orthomclBlastParser};
	system "$orthomclBlastParser $outDir/blast_all-v-all_output.blastn $outDir/compliantFasta >> $outDir/similarSequences.txt";
}

#==============================================
sub step9 {
    print "In step 9\n";

	open( CONFIG9, ">$outDir/config_step9" ) or die "can't open file $outDir/config_step9";

	print CONFIG9 "dbVendor=mysql
dbConnectString=dbi:mysql:mysql_local_infile=1:o_" . $tablesNamesSuffix . ":" . $mysqlHostName . "
dbLogin=". $mysqlUserName. "
dbPassword=" . $mysqlPassword ."
similarSequencesTable=SimilarSequences" . $tablesNamesSuffix;

	my $orthomclLoadBlast = $Config->{orthomcl}->{orthomclLoadBlast};
	system "$orthomclLoadBlast $outDir/config_step9 $outDir/similarSequences.txt";

	close(CONFIG9);
}

#==============================================
sub step10 {
    print "In step 10\n";

	open( CONFIG10, ">$outDir/config_step10" ) or die "can't open file $outDir/config_step10";

	print CONFIG10 "dbVendor=mysql
dbConnectString=dbi:mysql:mysql_local_infile=1:o_" . $tablesNamesSuffix . ":" . $mysqlHostName . "
dbLogin=". $mysqlUserName. "
dbPassword=" . $mysqlPassword ."
similarSequencesTable=SimilarSequences" . $tablesNamesSuffix . "
orthologTable=Ortholog" . $tablesNamesSuffix . "
inParalogTable=InParalog" . $tablesNamesSuffix . "
coOrthologTable=CoOrtholog" . $tablesNamesSuffix . "
interTaxonMatchView=InterTaxonMatch" . $tablesNamesSuffix . "
percentMatchCutoff=" . $percentMatchCutoff . "
evalueExponentCutoff=-5";

	close(CONFIG10);
	#my $orthomclPairs = $Config->{orthomcl}->{orthomclPairs};
	#system "$orthomclPairs $outDir/config_step10 $outDir/orthomclPairs.log cleanup=yes";
	my $cmd = "perl -w $scriptsDir/orthomclPairs.pl $outDir/config_step10 $outDir/orthomclPairs.log cleanup=yes $mysqlLibPath";
	print "orthomclPairs command: $cmd\n";
	system $cmd;

}

#==============================================

sub step11{
    print "In step 11\n";

	system "perl -w $scriptsDir/orthomclDumpPairsFiles.pl $outDir/config_step10 $outDir $mysqlLibPath";
}

#==============================================

sub step12{
    print "In step 12\n";

	my $mcl = $Config->{orthomcl}->{mcl};
	system "$mcl $outDir/mclInput --abc -I $inflationParam -o $outDir/mclOutput -force-connected y";
}

#==============================================

sub step13{
    print "In step 13\n";

	my $orthomclMclToGroups = $Config->{orthomcl}->{orthomclMclToGroups};
	system "$orthomclMclToGroups $genusName 0 < $outDir/mclOutput > $outDir/groups.txt";
}