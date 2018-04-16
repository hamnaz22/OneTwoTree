use strict;
use Bio::SeqIO;

if (@ARGV < 3 ){
	die "usage:
	(1) input file (groups.txt - orthomcl's output)
	(2) output directory - will contain a file for each group
	(3) all sequences fasta file";
}

my ($input, $outputDir, $fastaAll) = @ARGV;
my ($line,$groupName, @group, $GI, %seqsHash,$fileName, $seqCoverage, $numberOfSeqInCluster,  $totalNumberOfSequences, $minSeqCoverage, $header, $sequence);
my $c = 0;
my ($genusName) = ($fastaAll =~ m/([^\/]+)-allseq/);
%seqsHash = ();

open(IN,"<$input") or die "can't open file $input";


fillSeqsHash();


unless (-d $outputDir){
	mkdir($outputDir);
}

my @lines = <IN>;

print "Genus name: $genusName\n";

foreach my $line (@lines){
	chomp($line);
	if ($line =~ m/.+: (.+)/ || $line =~ m/(.+)/){
		$groupName = $genusName."-".($c < 10 ? "0$c" : $c);
		$c++;
		$fileName = "$outputDir/$groupName.txt";
		@group = split(/ /, $1); #each group contain IDs in the form:  taxon_id|sequence_id
		print "group is $line \n";
		my %organismsHash = ();
		my ($header, $taxonID, $numOfOrganismsInCluster);
			
		#count how many organisms exist in the group (cluster)
		foreach my $ID (@group){
		    #print "in ID $ID\n";
			#if($ID =~ m/\d+\|(\d+)/ || $ID =~ m/(\d+)/){
			if($ID =~ m/\d+\|([A-Z{2}\d]+\.\d+)/ || $ID =~ m/([A-Z{2}\d]+\.\d+)/){
				$GI = $1;			
				$header = $seqsHash{$GI}{"header"};
				($taxonID) = ($header =~ m/taxonid\|(\d+)\|/);
				if(!defined $taxonID){
					print "GI: $GI ID: $ID\n";
				}
				$organismsHash{$taxonID} = "1";
			}
		}

		$numOfOrganismsInCluster = scalar(keys(%organismsHash));
        print "organisms in cluster = $numOfOrganismsInCluster\n";

		if ($numOfOrganismsInCluster > 3){ #take only large enough clusters
		    print "Found large enough cluster $numOfOrganismsInCluster\n";
			open(OUT, ">$fileName") or die "can't open file $fileName";
			foreach my $ID (@group){
				#if($ID =~ m/\d+\|(\d+)/ || $ID =~ m/(\d+)/){
				if($ID =~ m/\d+\|([A-Z{2}\d]+\.\d+)/ || $ID =~ m/([A-Z{2}\d]+\.\d+)/){
					$GI = $1;	
					$header = $seqsHash{$GI}{"header"};
					$sequence = $seqsHash{$GI}{"sequence"};
					print OUT "$header\n";
					print OUT "$sequence\n";
				}
			}
			close(OUT);
		} else {
		    print "Not enough organisms for cluster $numOfOrganismsInCluster\n";
		}
	}
}

close(IN);


#insert all the sequences to seqsHash
sub fillSeqsHash {
	my $in = Bio::SeqIO->new("-file"=>"<$fastaAll","-format"=>"Fasta");	#open all-seqs fasta file
	my ($line, $GI, $header, $sequence, $restHeader, $oldGI, $fastaSeqObj);
	while (defined($fastaSeqObj = $in->next_seq())) { #iterate over all the sequences
		$header = getHeader($fastaSeqObj); #get the header of the fasta entry
	
		#unique gi:		
		($GI) = ($header =~ m/gi\|(.+)\|taxonid/);
		#print "GI=$GI\n";
		
		# index the headers and the sequences in the hash by GI --> header & sequence
		$seqsHash{$GI}{"header"} = ">$header";
		$seqsHash{$GI}{"sequence"} = $fastaSeqObj->seq();
	}
}


sub getNumOfOrganisms{
	my %organismsHash = ();
	my $header;
	my @keys = keys(%seqsHash);
	foreach my $key (@keys){
		$header = $seqsHash{$key}{"GI"};
		my ($taxonID) = ($header =~ m/taxonid\|(\d+)\|/);
		$organismsHash{$taxonID} = "1";
	}
	return scalar(keys(%organismsHash));
}



sub getHeader{
	my $fastaSeqObj = $_[0];
	my $id = $fastaSeqObj->id();
	my $desc = $fastaSeqObj->desc();
	my $header = $id." ".$desc;
	return $header;
}