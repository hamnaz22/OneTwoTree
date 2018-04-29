use strict;
use Bio::SeqIO;

#inputs: (1) fasta file of sequences of multiple species of the same genus
#		 (2) directory for the output files - each file for the sequences of the same species
if(@ARGV != 2){
	die "Ivalid number of arguments @ARGV usage:
	(1) fasta file of sequences of multiple species of the same genus
	(2) directory for the output files - each file for the sequences of the same species";
}

my ($inFile, $dir) = @ARGV;
my $in = Bio::SeqIO->new("-file"=>"<$inFile","-format"=>"Fasta");	#open all-seqs fasta file

mkdir($dir) unless (-d $dir);

my %speciesHash = ();
my ($path, $fileName, $out, $fastaSeqObj);

while (defined($fastaSeqObj = $in->next_seq())){
	my $header = getHeader($fastaSeqObj);
	if ($header =~ m/gi\|.*\|taxonid\|(\d*)\|organism\|/){ #if header, read the species name
		$fileName = $1;
		$path = $dir.'/'.$fileName.".fasta";
		if (!exists $speciesHash{$fileName}){ #the taxon file doesn't exist
			$out = Bio::SeqIO->new("-file" => ">$path", "-format" => "Fasta");
			$speciesHash{$fileName}=1;
		}
		else{ #the taxon file exist
			$out = Bio::SeqIO->new("-file" => ">>$path", "-format" => "Fasta");
		}
		$out->write_seq($fastaSeqObj);
	}
}

#get header from sequence object
sub getHeader{
	my $fastaSeqObj = $_[0];
	my $id = $fastaSeqObj->id();
	my $desc = $fastaSeqObj->desc();
	my $header = $id." ".$desc;
	return $header;
}
