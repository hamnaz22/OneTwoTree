package Guidance;
use strict;
use Bio::SeqIO;
use Config::Tiny;

sub runGuidance{
	if(@_ < 4){
			die "usage:
			(1) sequences of core MSA (unaligned)
			(2) fragments (unaligned)
			(3) output directory
			(4) config .ini file
			";
	}
	my ($sequencesOfCoreMSA, $fragments, $outDir, $configFile) = @_;
	
	# Create a config
	my $Config = Config::Tiny->new;
	# Open the config
	$Config = Config::Tiny->read( $configFile );
	# Reading properties
	my $rootproperty = $Config->{_}->{rootproperty};
	my $guidance = $Config->{guidance}->{guidance};
	
	# With adjust direction
	print "Calling perl $guidance --seqFile $sequencesOfCoreMSA --msaProgram MAFFT --seqType nuc --outDir $outDir --MSA_Param \"\\-\\-adjustdirection \\-\\-addfragments $fragments \\-\\-multipair\"\n";
	system "perl $guidance --seqFile $sequencesOfCoreMSA --msaProgram MAFFT --seqType nuc --outDir $outDir --MSA_Param \"\\-\\-adjustdirection \\-\\-addfragments $fragments \\-\\-multipair\"";

}

# Filter out low-scored sequences from the MSA of ITS:
# Get the output of guidance and filter out the sequences that received score less than 0.6 and less than mean-score minus 3*standard-deviation
#NOTE: this method is not optimal, better methods should be tested: http://en.wikipedia.org/wiki/Chauvenet%27s_criterion
sub checkGuidance{
	if(@_ < 3){
		die "usage:
		(1) input MSA
		(2) output MSA (to contain the higher-scored sequences
		(3) the output directory of guidance
		";
	}
	my ($inMSA, $outMSA, $guidanceDir) = @_;

	my $SD_threshold = 3;
	my @scores              = ();
	my %giScore              = ();
	my $str             = "";
	my $scorePerSeqFile = "$guidanceDir/MSA.MAFFT.Guidance_res_pair_seq.scr_with_Names";

	my $in = Bio::SeqIO->new("-file"=>"<$inMSA","-format"=>"Fasta");
	my $out = Bio::SeqIO->new("-file"=>">$outMSA","-format"=>"Fasta");
	my %seqHash = ();

	# index the input sequences in a hash of gi-->seqObj
	while (defined(my $seqObj = $in->next_seq())) {
		my $id = $seqObj->id(); 
		my ($gi) = ($id =~ m/gi\|(\d+)\|/);
		$seqHash{$gi} = $seqObj;	
	}

	open( SCR, $scorePerSeqFile ) or die "can't open file $scorePerSeqFile";

	while ( my $line = <SCR> ) {
		chomp($line);
		if ( $line =~ m/gi\|(\d+)\|.*\t([01]\.\d+)$/ ) {
			my $gi    = $1;
			my $score = $2;
			push( @scores, $score );
			$giScore{$gi} = $score;
		}
	}

	my $mean = mean(@scores);
	my $sd = standard_deviation(@scores);
		
	foreach my $gi (keys(%giScore)){
		if ($giScore{$gi} <= 0.6 && $giScore{$gi} < $mean - $SD_threshold * $sd){
			$str .= "$gi - $giScore{$gi}\n";		
		}
		else{
			if (! exists $seqHash{$gi}){
				print "!!! $gi\n";
			}
			else{
				$out->write_seq($seqHash{$gi});
			}
		}
	}
	print "CheckGuidance: ".scalar(@scores)." sequences, average: $mean, standard deviation: $sd. Filtered out:
	$str\n";

	close(SCR);


	sub mean{
		my (@numbers) = @_;
		
		my $sum = 0;
		foreach my $n (@numbers){
			$sum += $n;	
		}
		return (@numbers > 0 ? $sum / scalar(@numbers) : 0);
	}

	sub standard_deviation {
		my (@numbers) = @_;

		#Prevent division by 0 error in case you get junk data
		return undef unless ( scalar(@numbers) );

		# Step 1, find the mean of the numbers
		my $total1 = 0;
		foreach my $num (@numbers) {
			$total1 += $num;
		}
		my $mean1 = $total1 / ( scalar @numbers );

		# Step 2, find the mean of the squares of the differences
		# between each number and the mean
		my $total2 = 0;
		foreach my $num (@numbers) {
			$total2 += ( $mean1 - $num )**2;
		}
		my $mean2 = $total2 / ( scalar @numbers );

		# Step 3, standard deviation is the square root of the
		# above mean
		my $std_dev = sqrt($mean2);
		return $std_dev;
	}
}



1;
