#!/bin/bash/perl -w
use strict;

#Input: (1) in fasta file name (2) out nexus file name

# SETUP
die "Usage: fasta2nexus.pl <in.fa.file MSA> <out.nxs.file>" unless @ARGV == 2;

# FIELD
my $fastaFile = $ARGV[0];
my $nexusFile = $ARGV[1];
my $taxNum = 0;
my $charNum = -1;
my %seq = ();

#print"line 16 $fastaFile\n$nexusFile\n";
# RUN
&readFasta($fastaFile);
&writeNexus($nexusFile);

# SUBROUTINES
sub readFasta{
	my $fileName = $_[0];
	open FASTA, $fileName or die "unable to open $fileName: $!";
	my ($seq_id, $seq_seq) = ('', '');
	while(<FASTA>){
		chomp;
		if(/>/){
			$seq{ $seq_id } = $seq_seq;
			($seq_id, $seq_seq) = ('', '');
			$seq_id = substr((split /\@/, $_)[0], 1);
			$taxNum++;
		}else{
			$seq_seq .= $_;
		}
	}
	close FASTA;
	$seq{ $seq_id } = $seq_seq;
	$charNum = length( $seq_seq );
	delete( $seq{''} );
}

sub writeNexus {
	my $outfile = $_[0];
	open NEXUS, "> $outfile" or die "unable to open $outfile $!";
	#print header
	print NEXUS "#NEXUS\n";
	print NEXUS "begin data;\n";
	print NEXUS "dimensions ntax=$taxNum nchar=$charNum;\n";
	print NEXUS "format datatype=DNA interleave=no gap=- missing=?;\n";
	print NEXUS "matrix\n";
	foreach(keys %seq){
		print NEXUS "$_\t";
		print NEXUS $seq{ $_ }."\n";
	}
	print NEXUS ";\n";
	print NEXUS "end;\n";
	close NEXUS;
}
