use strict;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;

if(@ARGV < 2){
	die "usage:
	(1) all-seq original fasta file (non-split)
	(2) output file - to contain all the sequences contain ITS";
}

my ($fastaAll, $outITSseqs) = @ARGV;
my ($fastaSeqObj, $header, $desc);

my $in = Bio::SeqIO->new("-file"=>"<$fastaAll","-format"=>"Fasta");	#open all-seqs fasta file
my $out = Bio::SeqIO->new("-file" => ">$outITSseqs", "-format" => "Fasta");

while (defined($fastaSeqObj = $in->next_seq())) { #iterate over all the sequences
	$header = getHeader($fastaSeqObj);
	($desc) = ($header =~ m/\|description\|(.*)$/);
	if($desc =~ m/ITS/ || $desc =~ m/internal transcribed spacer/ ){
		$out->write_seq($fastaSeqObj);
	}
}



sub getHeader{
	my $fastaSeqObj = $_[0];
	my $id = $fastaSeqObj->id();
	my $desc = $fastaSeqObj->desc();
	my $header = $id." ".$desc;
	return $header;
}