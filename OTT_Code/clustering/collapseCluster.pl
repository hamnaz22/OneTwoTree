use strict;
use Bio::SeqIO;
use Bio::DB::GenBank;
use Bio::Index::GenBank;

if ( @ARGV < 4) {
	die "usage: 
	(1) directory contains onlyDesc files
	(2) GenBank file of this genus
	(3) directory contains the sequence file of each cluster
	(4) all sequences fasta file
	";
}

my ( $onlyDescDir, $GBFile, $seqsDir, $allSeqFastaFile ) = @ARGV;

my %originalToUniqeGI;

opendir( DIR, "$onlyDescDir" ) or die "can't open $onlyDescDir";
my @onlyDescFiles = grep { $_ ne '.' && $_ ne '..' } readdir(DIR);

# for each "onlyDesc" file in the directory, collapse it
foreach my $file (@onlyDescFiles) {
	if ( $file !~ m/clp/ && $file ne "appendedClusters.txt" ) {
		my $inPath = "$onlyDescDir/$file";
		my $clpFileName;
		if ( $file =~ m/(\d+)/ ) {
			$clpFileName = "clp$1" . ".txt";
		}
		else {
			$clpFileName = "clp$file";
		}
		collapse( $inPath, $clpFileName, $GBFile );
	}
}

#take a single cluster and collapse its sequences by unifing duplicate feature names
sub collapse {

	my ( $cluster, $clpFileName, $GBFile) = @_;
	my ( %hash, $GI, $desc );
	%hash = ();
	open( IN, "<$cluster" ) or die "can't open file $cluster";
	my $gbIndex = Bio::Index::GenBank->new(-filename   => $GBFile . ".idx",	-write_flag => 1);
	$gbIndex->make_index($GBFile);

	while ( my $line = <IN> ) {
		
		$desc = getFeatures($line, $gbIndex);
		
		# hash - key: description, value: number of occurences
		if ( !exists $hash{$desc} ) {
			$hash{$desc} = 1;
		}
		else {
			$hash{$desc} = $hash{$desc} + 1;
		}
	}
	close(IN);

	my @descs = keys(%hash);
	my @sorted = sort { $hash{$b} <=> $hash{$a} } @descs;
	
	# ------------------------- write the clusters: -----------------------
	my $clpFilePath = "$onlyDescDir/$clpFileName";
	writeCollapsedCluster([@sorted], {%hash}, $clpFilePath);	
}

# get the description of the feature(s) of a sequence, according to the "feature" field in a split sequence, or the GenBank file in a non-split sequence
sub getFeatures{
	my ($line, $gbIndex) = @_;
	my $desc;
	# support various formats of header lines
	if ( $line =~ m/\|feature\|.*:(.*)\|old_gi/ ) {
		$desc = $1;
	}
	elsif ( $line =~ m/\|feature\|.*:(.*)/ ) {
		$desc = $1;
	}
	elsif ( $line =~ m/>(\d+)\|\d+\|(.*)/ ) {
		my $GI = $1;
		$desc = getFeaturesInfo( $GI, $gbIndex );
	}

	$desc = formatKeyWords($desc);
	return $desc;
}


#write a collapsed cluster to a file
sub writeCollapsedCluster{
	if(@_ < 3){
		die "usage:
		(1) array of descriptions
		(2) hash of description-->number of occurences
		(3) output path
		";
	}
	
	my ($arrRef, $hashRef, $clpFile) = @_;
	my @descs = @{$arrRef};
	my %hash = %{$hashRef};
	open( OUT, ">$clpFile" ) or die "can't open file $clpFile";
	foreach my $desc (@descs) {
		my $n = $hash{$desc};
		if ( length($n) < 2 ) {
			print OUT "> " . $n . "\t\tX\t$desc\n";
		}
		else {
			print OUT "> " . $n . "\tX\t$desc\n";
		}
	}
	close(OUT);
}


# format the free-text to get a better collapse
sub formatKeyWords {
	my ($desc) = @_;

	$desc =~ s/internal transcribed spacer 1/ITS1/ig;
	$desc =~ s/ITS 1/ITS1/ig;
	$desc =~ s/internal transcribed spacer 2/ITS2/ig;
	$desc =~ s/ITS 2/ITS2/ig;
	$desc =~ s/ribosomal RNA/rRNA/ig;
	$desc =~ s/trnH-psbA/psbA-trnH/ig;
	$desc =~ s/tRNA-leu/trnL/ig;
	$desc =~ s/intergenic spacer region/intergenic spacer/ig;
	$desc =~ s/rpLl16/rpL16/ig;
	$desc =~ s/genomic DNA \| //ig;
	$desc =~ s/(.+)\s* \(\1\)/$1/ig;
	$desc =~ s/coxi/cox1/ig;
	$desc =~ s/, and/ and/ig;
	
	#delete duplicates inside the feature (delimited by ';', ':' or ',')
	while ( $desc =~ s/(\|?\s*)(\S.+\S)\s*(;|:|,)\s*\2(\s*(\|?|$))/$1$2$4/ig ) {}

	#delete duplicate features:
	while ( $desc =~ s/((\S.+\S)\s*)\|\s*\2\s*(\||$)/$1$3/ig ) { }

	#clean white spaces in the start and end of the string
	$desc =~ s/^\s*(\S.+\S)\s*$/$1/i;
	
		$desc =~ tr/A-Z/a-z/;
	$desc =~ s/its/ITS/g;
	$desc =~ s/rrna/rRNA/g;
	return $desc;
}

# get the value of the sub-feature indicates the type of the feature object (e.g. which gene is it)
sub getFeatureValue{
	my ($featureObject) = @_;
	my $featureValue;
	for my $tag ( $featureObject->get_all_tags ) {
		if ( $tag eq "gene" ) {
			$featureValue = join( ' ', $featureObject->get_tag_values($tag) );
			return $featureValue;
		}
	}

	for my $tag ( $featureObject->get_all_tags ) {
		if ( $tag eq "product" ) {
			$featureValue = join( ' ', $featureObject->get_tag_values($tag) );
			return $featureValue;
		}
	}

	for my $tag ( $featureObject->get_all_tags ) {
		if ( $tag eq "note" ) {
			$featureValue = join( ' ', $featureObject->get_tag_values($tag) );
			return $featureValue;
		}
	}

	for my $tag ( $featureObject->get_all_tags ) {
		if ( $tag eq "mol_type" ) {
			$featureValue = join( ' ', $featureObject->get_tag_values($tag) );
			return $featureValue;
		}
	}
 }

# get the features of a sequence according to the GenBank file, and concatenate all of them
sub getFeaturesInfo {
	my ( $GI, $gbIndex ) = @_;
	my ( @features, $featuresInfo, $featureValue, $featureType );
	my $gbSeqObj = $gbIndex->fetch($GI);    # get the entry from GenBank
	@features =	  $gbSeqObj->get_SeqFeatures;    #get the sequence's features (cds, source, gene....)
	
	foreach my $featObj (@features) {
		$featureValue = getFeatureValue($featObj);
		$featureType  = $featObj->primary_tag;

		#print "$featureType: $featureValue\n";
		if ( !defined $featuresInfo ) {
			$featuresInfo = $featureValue;
		}
		else {
			$featuresInfo .= " | " . $featureValue;
		}
	}

	return $featuresInfo;
}

