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
";
}

my ($fastaAllSequences, $tablesNamesSuffix, $percentMatchCutoff, $inflationParam, $scriptsDir, $gbFile, $outDir, $configFile, $filterBlastRatio) = @ARGV;    #default: $percentMatchCutoff = 50, $inflationParam = 1.5
my (@seqsByOrganismFiles, $taxonCode, $numOfSequnces, $genusName);


# Create a config
my $Config = Config::Tiny->new;
# Open the config
$Config = Config::Tiny->read( $configFile );
# Reading properties
my $rootproperty = $Config->{_}->{rootproperty};
my $mysqlHostName = $Config->{ploidb_mysql}->{hostname};
my $mysqlUserName = $Config->{ploidb_mysql}->{username};
my $mysqlPassword = $Config->{ploidb_mysql}->{password};


#run_blastClust();


# parse the output:
my $onlyDesc = $outDir."_onlyDesc";
my $seqs = $outDir."_seqs";

system "perl -w $scriptsDir/getSequencesFromClusters.pl $outDir/BC_groups.txt $seqs $fastaAllSequences"; #get the sequences

system "perl -w $scriptsDir/getOnlyDesc.pl $seqs $onlyDesc"; #get descriptions of each cluster

system "perl -w $scriptsDir/collapseCluster.pl $onlyDesc $gbFile $seqs $fastaAllSequences n"; #collapse the clusters

system "perl -w $scriptsDir/appendClusters.pl $onlyDesc"; #append the collapsed clusters to one file - easier to read
