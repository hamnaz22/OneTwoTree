package OneSeq;

use strict;
use Bio::SeqIO;

####### Split2Species ###################################################################### 
# Recieves directories containing Orthomcl clusters and separates them according to species.
# Each sub-directory will contain fasta files of the separated species.
# Input: directory with sub-directories which contain sequences files
############################################################################################

sub Split2Species {
	my ($file,$dir)=@_;
	my $in = Bio::SeqIO->new("-file"=>"$file","-format"=>"Fasta");	#open all-seqs fasta file
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
}

#get header from sequence object
sub getHeader{
	my $fastaSeqObj = $_[0];
	my $id = $fastaSeqObj->id();
	my $desc = $fastaSeqObj->desc();
	my $header = $id." ".$desc;
	return $header;
}

sub BuildBlastHash {
	my ($blast,%combinationScores,$min,$i,$de_combinationScores,%giHash);
	my ($blastFile,$cluster)=@_;
	open (IN,"<$cluster") or die "Can't open $cluster";
	my @clusterLines=<IN>;
	my $gi=CreateGIArray([@clusterLines]);
	my @gi=@{$gi};	
	close (IN);
	%combinationScores=();
	%giHash=();
	$min=1000;
	
	foreach my $key (@gi){
		$giHash{$key}=1;
	}
	
	open (BLAST,"<$blastFile") or die "Can't open $blastFile";
	$blast=<BLAST>;
	while (defined $blast){
		#if (($blast=~m/\d+\|(\d+)\s\d+\|(\d+).+\s(\d+\.?\d?)$/)||($blast=~m/gi\|(\d+)\|taxonid\|\d+.+gi\|(\d+)\|taxonid\|\d+.+\s(\d+\.?\d?)$/)){ #adjust to fit accession number istead of GI number:
		if (($blast=~m/([A-Z{2}\d]+\.\d+)\|(\d+)\s([A-Z{2}\d]+\.\d+)\|(\d+).+\s(\d+\.?\d?)$/)||($blast=~m/gi\|([A-Z{2}\d]+\.\d+)\|taxonid\|\d+.+gi\|([A-Z{2}\d]+\.\d+)\|taxonid\|\d+.+\s(\d+\.?\d?)$/)){
			if (((exists ($giHash{$1})) || (exists($giHash{$2})))&& ($1!=$2)){
				($min,$de_combinationScores)=subHashBuilder($1,$2,$3,{%combinationScores},$min);
				%combinationScores=%{$de_combinationScores};
				
			}
		}
		$blast=<BLAST>;
	}
	close (BLAST);
	return ({%combinationScores},$min);
}

sub subHashBuilder {
	my ($a,$b,$c,$de_combinationScores,$min)=@_;
	my %combinationScores=%{$de_combinationScores};
	if (exists ($combinationScores{$a}{$b})){ #if there is more than one match of the combination-->sum it up
		$combinationScores{$a}{$b}=$combinationScores{$a}{$b}+$c;
		if ($combinationScores{$a}{$b}<$min){
			$min=$combinationScores{$a}{$b};
		}						
	}
	else {
		$combinationScores{$a}{$b}=$c;
		if ($combinationScores{$a}{$b}<$min){
			$min=$combinationScores{$a}{$b};
		}						
	}
	return ($min,{%combinationScores});	
}

sub CreateGIArray {
	my ($de_fastaLines)=@_;
	#print("DEBUG## de_fastaLines: $de_fastaLines \n");
	my @fastaLines=@{$de_fastaLines};
	my @gi=();
	foreach my $line (@fastaLines){
		#if ($line=~m/>gi\|(\d*)/){ Changed to recognize the accession id in the gi feild
		#print("DEBUG## line in fastaLines: $line \n");
		if ($line=~ m/gi\|([A-Z{2}\d]+\.\d+)\|/){   #($header =~ m/gi\|([A-Z\d]+).(\d+)\|/)
			push(@gi,$1);
			print("DEBUG## Push gi: $1 \n");
		}
	}
	return ([@gi]);
}

####### ReturnMaxBitScore ###################################################################### 
# Recieves references to a fasta file content, the combination hash and the genus name (or sub directory). 
# This subroutine calculates out of a multiple accession which is the best sequence, according to a normalized
# bit score (calculated using it's average).
# Input: directory containing blast all-v-all file.
############################################################################################

sub Exists{
	my ($subject,$query,$de_combinationScores)=@_;
	my %combinationScores=%{$de_combinationScores};
	return (exists ($combinationScores{$subject}{$query}));
}

sub BuildScores{
	my ($gi,$score,$de_scores)=@_;
	my %scores=%{$de_scores};
	if ((! exists ($scores{$gi}[0]))&& (! exists ($scores{$gi}[1]))){
		$scores{$gi}[0]=$score;
		$scores{$gi}[1]=1;	
	}
	else {
		$scores{$gi}[0]=$scores{$gi}[0]+$score;		
		$scores{$gi}[1]=$scores{$gi}[1]+1;
	}
	return({%scores});
}

#$gi: a vector of gi's of one species that exist in one cluster
#$de_combinationScores: hash reference for all-vs-all blast scores for all sequences in a cluster
#min = min score to assign for pairs pf sequences that do not appear in the blast hash
#$cluster: file name including the sequences in the cluster

sub ReturnMaxBitScore {
	my ($gi1,$de_combinationScores,$min,$cluster,$mode,$gi2)=@_;
	my @gi1=@{$gi1};	
	my @gi2;
	if (not defined $mode){
		$mode='first';
	}
	if (not defined $gi2){
		@gi2=@{$gi1};	
	}
	else {		
		@gi2=@{$gi2};	
	}
	my %combinationScores=%{$de_combinationScores};
	my ($i,$j,%scores);
	%scores=();
	my $de_scores={%scores};

	foreach (@gi1) {
	  print "##DEBUG gi1 array elements: $_\n";
	  print "$_\n";
	}

	foreach (@gi2) {
	  print "##DEBUG gi2 array elements: $_\n";
	  print "$_\n";
	}

	for ($i=0; $i<scalar(@gi1);$i++){
	
		print("##DEBUG i --> $i \n");
		
		for ($j=0; $j<scalar(@gi2);$j++){

			print("##DEBUG j --> $j \n");
			print("##DEBUG gi1[i] --> $gi1[$i] \n");
			
			#if ($gi1[$i]==$gi2[$j]){
			if ($gi1[$i] eq $gi2[$j]){	#changed since gi is no longer an integer but a string (accession number such as AZ2387263.1)
				print("##DEBUG $gi1[$i] == $gi2[$j] \n");
				next;
			}
				
			
			my $ij_exists = Exists($gi1[$i],$gi2[$j],{%combinationScores});
			print("##DEBUG ij_exists --> $ij_exists \n");
			#	my $ji_exists = Exists($gi2[$j],$gi1[$i],{%combinationScores});
			my $ij_score = $combinationScores{$gi1[$i]}{$gi2[$j]};
			#	my $ji_score = $combinationScores{$gi2[$j]}{$gi1[$i]};
			print("##DEBUG ij_score --> $ij_score \n");
						
			if ($ij_exists){
				print("##DEBUG ij_exists --> $ij_exists , score $ij_score\n");
				$de_scores=BuildScores($gi1[$i],$ij_score,$de_scores);
			}
			else {
				$de_scores=BuildScores($gi1[$i],$min,$de_scores);			
			}
		}
	}	

	print("DEBUG## Average using de_scores: $de_scores\n");
	
	%scores=%{$de_scores};
	my ($selectedGI,$max,$k);
	
	$max=Average($scores{$gi1[0]}[0],$scores{$gi1[0]}[1]);
	$selectedGI=$gi1[0];
	my $maxGI;
	
	for ($k=1; $k<scalar(@gi1); $k++){
		my $curScore = Average($scores{$gi1[$k]}[0],$scores{$gi1[$k]}[1]);
		if ($max<=$curScore){
			$maxGI=$selectedGI;
			if (($curScore==$max)&&($mode eq 'first')){ #if all are equal in distance from each other
				open (IN,"<$cluster") or die "Can't open $cluster";
				my @old_gi=($maxGI,$gi1[$k]);
				my @clusterLines=<IN>;
				close (IN);
				$gi1=CreateGIArray([@clusterLines]);
				$selectedGI=ReturnMaxBitScore([@old_gi],{%combinationScores},$min,$cluster,'second',$gi1);
			}
			$max=$curScore;
			$selectedGI=$gi1[$k];
			$mode='first';
		}
	}
	return ($selectedGI);
}

sub Average {
	my ($sum,$count)=@_;
	my $avg=$sum/$count;
	return ($avg)
}

1;