#####################################################################################################################################
#
#	Script : Create Nexus files for concatenated alignment
#	Description : The script creates the Nexus files required by MrBayes, for all input clusters and concat alignment
#
#	Input : cluster alignment fasta files, outgroups files, concat alignment fasta, concat report
#	Output : For each cluster and concat - nexus sequence file and nexus configuration file (later to be executed in MrBayes )
#
#	Parameters :
#	1. cluster list.txt (format explained below) 2. concat alignment.fasta 3. concat report.txt 4. concat outgroup.fasta ('no' if no outgroup available)
#	5. parameters string (see below) 6. concat configuration output.nex 7. concat sequence output.nex 8. Mrbayes output (no extension) 9. simplify taxa names mode (yes/no)
#
#####################################################################################################################################
#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;

#print("@ARGV\n\n");
print(scalar(@ARGV));
print("\n\n");
# Input

my $clusterList = $ARGV[0]; # a text file with the list of clusters fasta files and their outgroup fasta files
my $concatFile   = $ARGV[1];    # the concat fasta file
my $concatReport = $ARGV[2];    # the concat report text file
my $concatOutgroupFile = $ARGV[3];	# outgroup sequence\s file for concat
my $paramString = $ARGV[4];	# MrBayes parameters, see below. use give to use defaults
my $concatConfOut = $ARGV[5];	# concat config.nex
my $concatSeqOut = $ARGV[6];	# concat seq.nex
my $concatMbOut = $ARGV[7];		# concat MrBayes output
my $concatSimplify = $ARGV[8];	# yes or no, to simplify concat seq names
#Michal - enable more user control on MB params:
my $nchains_val = $ARGV[9];
my $samplefreq_val = $ARGV[10];
my $burninFrac_val = $ARGV[11];
my $checkFreq_val = $ARGV[12];
#New Clock parametrs:
my $clock_model = $ARGV[13]; # Unconstrained, Strict clock, Relaxed clock
# Unconstrained : remove the line prset brlenspr=....
# Strict clock: prset brlenspr=clock:uniform/birthdeath/coalescence
# Relaxed clock: 	prset brlenspr=clock:uniform;
#					prset clockvarpr=igr;  /Cpp/TK02/Igr
#my $branch_length_model = $ARGV[14];
my $branch_length_model = $ARGV[14];
print("@ARGV\n");


print("$samplefreq_val\n");

my $usage = "create_nexus_for_concat.pl
	1. cluster list
	2. concat alignment.fasta 
	3. concat report  
	4. concat outgroup.fasta ('no' if no outgroup available)
	5. parameters string (put in \"\")
	6. concat configuration output.nex
	7. concat sequence output.nex
	8. Mrbayes output (no extension)
	9. simplify taxa names mode (yes/no)
	10. nchains_val
	11. samplefreq_val
	12. burninFrac_val
	13. checkFreq_val
	14. clock_model
	15. branch_length_model
\n";

print("1\n");

#unless ( scalar(@ARGV) == 15  || scalar(@ARGV) == 19) {
#	die $usage;
#}

print("2\n");

### get parameters input
my %paramHash = (	
	"ngen" => 2000000,
	"nruns" => 2,
	#"nchains" => 4,
	"nchains" => $nchains_val,
	#"samplefreq" => 2000,
	"samplefreq" => $samplefreq_val,
	"diagnfreq" => 10000,
	"mcmcdiagn" => "Yes",
	"relburnin" => "yes",
	#"burninfrac" => 0.25,
	"burninfrac" => $burninFrac_val,
	#"Stoprule" => "Yes",
	#"Stopval" => 0.01,
	"conformat" => "simple",
	"Checkpoint" => "Yes",
	#"Checkfreq" => 10000
	"Checkfreq" => $checkFreq_val
	);

print("3\n");
# param string should include one or more of the parameters, separated by spaces
# e.g : "ncahins=6 ngen=2000000"
my @params = split(/ /,$paramString);
foreach my $p (@params){
	my ($paramName,$paramValue) = split(/=/,$p) or die "Bad parameter format for $p
	format should be: 'parameterName=parameterValue'\n";
	if (exists $paramHash{$paramName}){
		$paramHash{$paramName} = $paramValue;
	}
	else{
		die "Unknown parameter '$paramName'\n";
	}
}

my $mcmcString = "";
foreach my $p ("nruns","nchains","ngen","samplefreq","diagnfreq","mcmcdiagn","relburnin","burninfrac","Checkpoint","Checkfreq"){
	$mcmcString .= $p."=".$paramHash{$p}." ";
}

my $sumtString = "";
foreach my $p ("relburnin","burninfrac","conformat"){
	$sumtString .= $p."=".$paramHash{$p}." ";
}


##### Create nexus for clusters

### Read input cluster list
open( LIST, "< $clusterList" );
open(my $listLines, "<", $clusterList)
  or die "Can't open cluster list file $clusterList";
my @listLines = <LIST>;


### Create Nexus for list clusters

my @outgroupFilesList; # for later use in creating concat nexus
my $listClusterCount = 0; # counts lines for later control
print " CHCHCHCHC ListClusterCount variable  = $listClusterCount.......$clusterList ";
foreach my $listLine (@listLines) {
	print "$listLine  CMXMCMXCMXCMXCMXMC\n";
	chomp $listLine;
	$listClusterCount++;
	print " CHCHCHCHC ListClusterCount variable  = $listClusterCount";
# list line format should be: cluster.fasta;outgroup.fasta ('no' if not available);config.nex output;seq.nex output;mrbayes output;simplify(yes\no)
	if ( $listLine !~ m/([^;]+);([^;]+);([^;]+);([^;]+);([^;]+);(yes|no)/ ) {
		print
"WARNING : line $listLine does not match the expected line format. skipping line\n";
	}
	else {    # if line is OK
		my $msaFile      = $1;
		my $outgroupFile = $2;
		push (@outgroupFilesList, $outgroupFile);
		my $outConf      = $3;
		my $outSeq       = $4;
		my $mrBayesOut   = $5;
		my $simplify     = $6;
		create_nexus_for_cluster( $msaFile, $outgroupFile, $outConf, $outSeq,
			$mrBayesOut, $simplify );
	}
}
close (LIST);

##### Create nexus for concat

my $concatToConvert; # fasta file that will be converted to nexus

### Simplify seq names in concat file and report file, if option is enabled

my $simpleReport = $concatReport;
if ($concatSimplify eq "yes"){
	# simplify concat
	my $simpleConcat = $concatFile;
	$simpleConcat =~ s/(.+)\.fasta/$1_simple\.fasta/; 
	$concatToConvert = simplifyFastaHead ($concatFile, $simpleConcat);
		
	# simplify report
	$simpleReport =~ s/(.+)\.txt/$1_simple\.txt/;
	open (REP1, "< $concatReport") or die "Can't open report file";	# input report
	open (REP2, "> $simpleReport") or die "Can't open simple report file"; # simple report
	my @reportLines = <REP1>;
	foreach my $reportLine (@reportLines){
		if ($reportLine =~ m/^.*gi\|(\d+)\|taxonid\|\d+\|organism\|([^\s]+) ([^\|]+)[^\t]+(\t.+)/){
			print REP2 "$1"."_$2"."_$3$4\n";
		}
		else{
			print REP2 "$reportLine";
		}
	}
	close (REP1); close (REP2);
	
}
else{
	$concatToConvert = $concatFile;	
}

### Create concat seq.nex

#system( 'perl', "/groups/itay_mayrose/shiranabad/Maya/ploiDB/SVN/mrbayes/convert_fasta2nexus.pl",	# CHANGE PATH!!!
system( 'perl', "/bioseq/oneTwoTree/ott_scripts/convert_fasta2nexus.pl",	# CHANGE PATH!!!
	$concatToConvert, $concatSeqOut ) == 0
	or die "Can't convert $concatToConvert to Nexus !";

### Create outgroup block
my $outgroupBlock;
if ($concatOutgroupFile ne "no"){
	$outgroupBlock = "constraint ingroup partial = ";
	
	# get outgroup headers
	my @outgroupHeads; # stores outgroup names
	open (COG, "< $concatOutgroupFile") or die "Can't open concat outgroup file $concatOutgroupFile";
	my @outgroupLines = <COG>;
	foreach my $outgroupLine (@outgroupLines){
		if ($outgroupLine =~ m/^>(.+)/){
			my $outgroupHead = $1;
			if ($concatSimplify eq "yes"){
				unless ($outgroupHead =~ m/^>.*gi\|(\d+)\|taxonid\|\d+\|organism\|([^\s]+) ([^\|]+)/){
					warn "Could not simplify outgroup name $outgroupHead\n"
				}
				push (@outgroupHeads, "$1" . "_$2" . "_$3");	
			}
			elsif ($concatSimplify eq "no"){
				push (@outgroupHeads, $outgroupHead);
			}
			
		}
	}
	unless ( @outgroupHeads > 0 ) {
		die "Outgroup not specified for concat. Could not get outgroup from file $concatOutgroupFile";
	}
	
	# insert ingroup taxa to outgroup block
	open (SEQ, "< $concatSeqOut") or die "Can't open $concatSeqOut";
	my @nexSeqLines = <SEQ>;
	foreach my $nexSeqLine (@nexSeqLines) {
		if ( $nexSeqLine =~ m/.+\t/ ) {    # check if line is a sequence line
			my $outgroup = 0;    # 0 for heads of ingroup, 1 for GIs of outgroup
			foreach my $head (@outgroupHeads)
			{   # check agains all outgroup headers
				if ( $nexSeqLine =~ m/$head/ ) {
					$outgroup = 1;
					last;
				}
			}
			if ( $outgroup == 0 )
			{   # after checking against all outgroup headers
				$nexSeqLine =~ m/(.+)\t/;
				$outgroupBlock .= "$1 ";
			}
		}
	}
	
	
	# insert outgroup taxa to outgroup block
	my $outgroupJoin = join( " ", @outgroupHeads );
	$outgroupBlock .= ": $outgroupJoin;\n";
}

close(SEQ);

### Create partition block

my $partitionBlock; 
my %partitions;
# collect all partitions and store in hash
open (REP, "< $simpleReport") or die "Can't open report file $simpleReport";
my @reportLines = <REP>;
foreach my $reportLine (@reportLines){
	if ($reportLine =~ m/[^\t]+\t(\d+)\t(\d+)\t(\S+)/){	# if line is not header
		my $start = $1+1; # start position in concat alignment
		my $end = $2;	  # en position in concat alignment
		my $cluster = $3;
		if ($cluster ne "?"){		# if line contains cluster name
			$cluster .= "_$start\_$end";
			$cluster =~ s/.+\\(.+)\.fasta/$1/;
			unless (exists $partitions{$cluster}) {	 # if cluster name has not been adde before
				$partitions{$cluster} = [$start,$end];
			}
		}
	}
	
}

# write partitions block according to hash

my @clusterKeys = keys(%partitions);
my $hashClusterCount = scalar(@clusterKeys);
foreach my $key (@clusterKeys){
	$partitionBlock .= "charset $key = "."$partitions{$key}[0]"."-"."$partitions{$key}[1]".";\n";
}

if ($hashClusterCount == $listClusterCount){
	$partitionBlock .= "partition clusters = $hashClusterCount : ".join(",",@clusterKeys).";\n"."set partition = clusters;\n";
}
else{
	die "Number of clusters in input list ($listClusterCount) does not match the number assumed from report file ($hashClusterCount)"
}

### create MrAIC block
my $mrAicBlock;
my $MrAicFile;
my $partitionIndex = 0; # used in the MrAIC block 
foreach my $key (@clusterKeys){
	$key =~ s/(.+)_\d+_\d+/$1/;
	$partitionIndex++;
	# find location of MrAIC output file
	open( LIST, "< $clusterList" ) or die "Can't open cluster list file $clusterList";
	my @listLines = <LIST>;
	my $path;
	foreach my $listLine (@listLines){
		if ($listLine =~ m/(.+$key);/){ # find the line that refers to the cluster		if ($listLine =~ m/$key;[^;]+;[^;]+;([^;]+)/){ # find the line that refers to the cluster
			$path = dirname($1); # path to the output nexus file, where the MrAIC output should also be
			#$path = "/groups/itay_mayrose/michaldrori/MSA_JmodelTree/output/Asparagus/concat/";
			last;
		}
	}
	close (LIST);
	
	# find the right file in the location
	opendir( DIR,  $path) or die "can't open $path for cluster $key"; 
	my @fileNames = grep { $_ ne '.' && $_ ne '..' } readdir(DIR);
	foreach my $file (@fileNames){
		#MD if ($file =~ m/MrAIC.txt$/){
		if ($file =~ m/mb_BLK.txt$/){
			$MrAicFile = "$path/$file";
			
		}
	}
	
	# parse MrAIC file to get the required block
	my $aiccBlock =""; # block for the current cluster
	open (MRAIC, "$MrAicFile") or die "Can't open MrAIC output file $MrAicFile for $key cluster";
	my $line      = <MRAIC>;
	while ( defined($line) ) {
		if ( $line =~ m/\BEGIN MRBAYES/ )
		{                      # find start line for the AICc block
			$line = <MRAIC>;
			while ( $line !~ m/END;/ ) {
				if ( $line =~ m/^ (Lset applyto=\(1\) .+)/ ) {
					$aiccBlock .= "$1\n";
					$line      = <MRAIC>;
				}
				if ( $line =~ m/^ (Prset applyto=\(1\) .+)/ ) {
					$aiccBlock .=  "$1\n";
					last;
				}
				$line = <MRAIC>;
			}
		}
		$line = <MRAIC>;
	}
	close (DIR);
	close (MRAIC);
	
	# Change the AIC block to the jModelTest Block
	$aiccBlock =~ s/applyto=\(1\)/applyto=($partitionIndex)/g; # insert partition index to block
	$mrAicBlock .= $aiccBlock; # add block for current cluster to general block
}


### Create concat config.nex
open (CONF, "> $concatConfOut") or die "Can't open configuration output file $concatConfOut";
print CONF "#NEXUS\n";
print CONF "begin mrbayes;\n";
print CONF "set autoclose=yes nowarn=yes;\n";
print CONF "execute $concatSeqOut;\n";
if ($concatOutgroupFile ne "no"){
	print CONF "$outgroupBlock";
}
print CONF "$partitionBlock";
print CONF "$mrAicBlock";
#MD#print CONF "prset brlenspr=clock:uniform;\n";
#prset brlenspr=clock:Strict;
#prset treeagepr=fixed(1);
#prset clockvarpr=igr;

if ($clock_model eq 'tk02' or $clock_model eq 'cpp' or $clock_model eq 'igr'){ # $branch_length_model
	print CONF "prset brlenspr=clock:$branch_length_model;\n";
	print CONF "prset clockvarpr=$clock_model;\n";
} elsif ($clock_model ne 'Unconstrained'){
	print CONF "prset brlenspr=clock:$clock_model;\n";
}

#prset brlenspr=clock:uniform;
#prset clockvarpr=igr;

if ($concatOutgroupFile ne "no"){
	print CONF "prset topologypr = constraints(ingroup);\n";
}

print CONF "mcmc $mcmcString file=$concatMbOut;\n";
print CONF "sumt $sumtString;\n";
print CONF "end;";


##### sub create nexus for cluster
sub create_nexus_for_cluster {
	my ( $msaFile, $outgroupFile, $outConf, $outSeq, $mrBayesOut, $simplify ) =
	  @_;

	# Create input file handles

# Create output file handles - Nexus files will be created in the MSA input location
	my $path        = dirname($outConf);
	my $clusterName = basename($msaFile);
	$clusterName =~ s/(.+)\.fasta/$1/;    # remove extension (.fasta)
	my $seqName = "$clusterName" . "_seq\.nex";
	print "$outConf\n";
	open( CONF, "> $outConf" )
	  or die "Can't open output nexus configuration file '$outConf'";    # configuration Nexus file

### Create sequence Nexus file

	my $fastaToConvert;    # the file to be later converted to Nexus format

# Simplify description of MSA fasta and outgroup file if 'simplify' option is enabled
	if ( $simplify eq "yes" ) {

		# MSA
		my $msaName = basename($msaFile);
		$msaName =~ m/([^\.]+)\..+/;
		my $simpleMSA = $path . "/$1" . "_simple.fasta";
		$fastaToConvert = simplifyFastaHead( $msaFile, $simpleMSA );

		# outgroup
		if ($outgroupFile ne "no"){
			my $outGroupName = basename($outgroupFile);
			$outGroupName =~ m/([^\.]+)\..+/;
			my $simpleOutgroup = $path . "/$1" . "_simple.fasta";
			$outgroupFile = simplifyFastaHead( $outgroupFile, $simpleOutgroup );
		}

	}
	elsif ( $simplify eq "no" ) {
		$fastaToConvert = $msaFile;
	}

	# Convert fasta to nexus using external script
	#system( 'perl', "/groups/itay_mayrose/shiranabad/Maya/ploiDB/SVN/mrbayes/convert_fasta2nexus.pl",
	system( 'perl', "/bioseq/oneTwoTree/ott_scripts/convert_fasta2nexus.pl",	
		$fastaToConvert, $outSeq ) == 0
	  or die "Can't convert $fastaToConvert to Nexus !";

### Get outgroup(s)
	my $outgroupBlock;
	if ($outgroupFile ne "no"){
		$outgroupBlock = "constraint ingroup partial = ";	# the actual text block to be added to the configuration file
		    	
		# get headers of outgroup taxon\taxa
		open( OUTGROUP, "< $outgroupFile" ) or die "Can't open outgroup input file $outgroupFile";
		my @outgroupHeads;    # stores headers from the outgroup file (one or more)
		my @outgroupLines = <OUTGROUP>;
		foreach my $outgroupLine (@outgroupLines) {
			chomp $outgroupLine;
			if ( $outgroupLine =~ m/^>(.+)/ ) {
				push( @outgroupHeads, $1 );
			}
		}
	
		unless ( @outgroupHeads > 0 ) {
			warn "Outgroup not specified in cluster $clusterName. Could not get outgroup from file $outgroupFile\n";
		}
	
		# create constraint ingroup:outgroup block
		# Ingroup
		open( SEQ, "< $outSeq" ) or die "Can't open Nexus sequence file";
		my @nexSeqLines = <SEQ>;
		foreach my $nexSeqLine (@nexSeqLines) {
			if ( $nexSeqLine =~ m/.+\t/ ) {    # check if line is a sequence line
				my $outgroup = 0;    # 0 for heads of ingroup, 1 for GIs of outgroup
				foreach my $head (@outgroupHeads)
				{                    # check agains all outgroup headers
					if ( $nexSeqLine =~ m/$head/ ) {
						$outgroup = 1;
						last;
					}
				}
				if ( $outgroup == 0 )
				{                    # after checking against all outgroup headers
					$nexSeqLine =~ m/(.+)\t/;
					$outgroupBlock .= "$1 ";
				}
			}
		}
		close(SEQ);
	
		# Outgroup
		my $outgroupJoin = join( " ", @outgroupHeads );
		$outgroupBlock .= ": $outgroupJoin;\n";
	}


### Run MrAIC


### Parse MrAIC output file - get AICc relevant NrBayes lines

	# Check if MrAIC output file was created
	#my $mraicOutput = "$msaFile\_codedNames_phy.phy\.MrAIC.txt";
	print("@@@@@@@@@@@@@@@@@@@@@@@@@".$msaFile."\n");
	$msaFile =~ s{\.[^.]*(?:\.fasta)?$}{};
	print("@@@@@@@@@@@@@@@@@@@@@@@@@".$msaFile."\n");
		my $mraicOutput = "$msaFile\_phy.phy.jmt_mb_BLK.txt";
	unless ( -e $mraicOutput ) {
		die "MrAIC output file not found !";
	}

	# Run on file lines and find the required ones
	open( MRAIC, "< $mraicOutput" );
	my $aiccBlock = "";        # will store the required MrBayes block (2 lines)
	my $line      = <MRAIC>;
	#$line      = <MRAIC>;
	print("START Line ".$line."\n");
	while ( defined($line) ) {
		#MD if ( $line =~ m/\[Mrbayes block for the best AICc model/ )
		if ( $line =~ m/\BEGIN MRBAYES;/ )
		{                      # find start line for the AICc block
			$line = <MRAIC>;
			while ( $line !~ m/END;/ ) {
				if ( $line =~ m/^ (Lset )applyto=\(1\) (.+)/ ) {
					$aiccBlock = "$aiccBlock" . "$1$2\n";
					$line      = <MRAIC>;
				}
				if ( $line =~ m/^ (Prset )applyto=\(1\) (.+)/ ) {
					$aiccBlock = "$aiccBlock" . "$1$2\n";
					last;
				}
				$line = <MRAIC>;
			}
		}
		$line = <MRAIC>;
		print($aiccBlock);
	}

	
	# check that the saved block has the expected formation
	if ( $aiccBlock !~ m/Lset.+\nPrset/ ) {
		die "Parsing failed ! bad AICc block formation";
	}

### Create configuration Nexus file

	# print to configuration Nexus file
	print CONF "#NEXUS  MDMDMDM\n";
	print CONF "begin mrbayes;\n";
	print CONF "set autoclose=yes nowarn=yes;\n";
	print CONF "execute $outSeq;\n";
	if ($outgroupFile ne "no"){
		print CONF "$outgroupBlock";
		print CONF "prset topologypr = constraints(ingroup);\n";
	}
	print CONF "$aiccBlock";
	if ($clock_model eq 'tk02' or $clock_model eq 'cpp' or $clock_model eq 'igr'){ # $branch_length_model
		print CONF "prset brlenspr=clock:$branch_length_model;\n";
		print CONF "prset clockvarpr=$clock_model;\n";
	} elsif ($clock_model ne 'Unconstrained'){
		print CONF "prset brlenspr=clock:$clock_model;\n";
	}
	print CONF "mcmc $mcmcString file=$mrBayesOut;\n";
	print CONF "sumt $sumtString;\n";
	print CONF "end;";

### Finish

	# close file handles
	close(OUTGROUP);
	close(SEQ);
	close(MSA1);
	close(MSA2);
	close(CONF);
	close(MRAIC);

	# delete files
	#system ("rm", "$fastaToConvert");
	#system ("rm", "$simpleOutgroup");

	# Subroutine : simplify fasta description
	sub simplifyFastaHead {
		my ( $inMSA, $outMSA ) = @_;
		open( MSA1, "< $inMSA" )  or die "Can;t open input MSA file";
		open( MSA2, "> $outMSA" ) or die "Can't create simplified MSA file";
		my @msaLines = <MSA1>;
		foreach my $msaLine (@msaLines) {
			if ( $msaLine =~
				m/^>?.*gi\|(\d+)\|taxonid\|\d+\|organism\|([^\s]+) ([^\|]+)/ )
			{
				my $simpleLine = ">$1" . "_$2" . "_$3\n";
				$simpleLine =~ tr/ '/_/
				  ; # remove spaces if there are any left (e.g in case of sub-species)
				print MSA2 "$simpleLine";
			}
			else {
				print MSA2 "$msaLine";
			}
		}
		close(MSA2);
		return ($outMSA);
	}
}

