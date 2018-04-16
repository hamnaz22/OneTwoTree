#given sequence fasta files of clusters, write only the descriptions of each file to the output directory
use strict;
use Bio::SeqIO;

if (@ARGV < 2 ){
	die "usage:
	(1) input directory - contain fasta files with the sequences of each cluster
	(2) output directory - will contain a file for each group
";
}

my ($inputDirClusterSequences, $outputDirOnlyDesc) = @ARGV;
my ($in, $header, $fastaSeqObj, $GI, $taxonID, $desc);

opendir( IN_DIR, $inputDirClusterSequences ) or die "can't open directory $inputDirClusterSequences";
my @SeqFiles = grep { $_ ne '.' && $_ ne '..' } readdir(IN_DIR);

mkdir($outputDirOnlyDesc) unless (-d $outputDirOnlyDesc);

foreach my $fileName (@SeqFiles){
	$in = Bio::SeqIO->new("-file"=>"<$inputDirClusterSequences/$fileName","-format"=>"Fasta");	#open all-seqs fasta file
	open(OUT, ">$outputDirOnlyDesc/$fileName") or die "can't open file $outputDirOnlyDesc/$fileName";
	
	while (defined($fastaSeqObj = $in->next_seq())) { #iterate over all the sequences
		$header = getHeader($fastaSeqObj); #get the header of the fasta entry

		($GI, $taxonID, $desc) = ($header =~ m/gi\|(.+)\|taxonid\|(\d+)\|.*\|description\|(.+)$/);
		print OUT ">$GI|$taxonID|$desc\n";
			
	}
	
	close(OUT);
}



sub getHeader{
	my $fastaSeqObj = $_[0];
	my $id = $fastaSeqObj->id();
	my $desc = $fastaSeqObj->desc();
	my $header = $id." ".$desc;
	return $header;
}