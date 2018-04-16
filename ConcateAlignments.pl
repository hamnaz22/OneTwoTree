use strict;
my $Aln_list=shift;
my $OutAln=shift;
my $Concate_report=shift;
my $PrintPerPos_Report=shift; #YES,NO (IF YES will print report: Pos_on_Concat_Aln,Species,Start_on_SeqBase_Aln,End_on_Base,Gene,Pos_On_Base_Aln,Reliaibility)
my $List_of_base_Aln_Info=shift; # OPTIONAL
my $Seq_codes_to_Names=shift; # OPTIONAL, to avoid give NA
my $Debug=shift; # YES will print QA prints, default=NO;
if (defined $Debug){$Debug=uc($Debug);}
else {$Debug="NO";}

my $Concate_report_Per_Pos=$Concate_report."_WithPerPos_Details";
open (LIST,$Aln_list) || die "Can't open list of alignment: '$Aln_list' $!";
my %Concate_Aln=();
my %Concate_Aln_Index=();
my $Concate_Aln_Hash_ref=\%Concate_Aln;
my $Concate_Aln_Index_Hash_ref=\%Concate_Aln_Index;
my $total_aln_length=0;
my %Seq_codes_to_Names=();
my @AlnLengths=(0); # each of the concatenated alignment will put here its length so we can creat blocks in the report for each one (including the missing ones)
my $AlnLengths_ref=\@AlnLengths;
if ($Seq_codes_to_Names ne "NA")
{
	open (SEQ_CODES,$Seq_codes_to_Names) || die "Can't open Seqs codes to names '$Seq_codes_to_Names' $!";
	while (my $line=<SEQ_CODES>)
	{
		chomp ($line);
		my ($seq_name,$code)=split(/\t/,$line);
		$Seq_codes_to_Names{$code}=$seq_name;
	}
}
while (my $aln_file=<LIST>)
{
	chomp ($aln_file);
	if (!-e $aln_file)
	{
		print "[WARNNING] $aln_file does not exist\n";
	}
	else
	{
		open (ALN,$aln_file) || die "Can't open alignment: '$aln_file' $!";
		my $gene_name=$aln_file;
		my @ans=validate_seq_names_unique($aln_file);
		if ($ans[0] ne "OK")
		{
			my @duplicates=@{$ans[1]};
			foreach my $dupl (@duplicates)
			{
				print "[ERROR] Sequence with name '$dupl' appear more than once on the file $aln_file\n";
			}
			die "ConcateAlignments.pl can't continue because the alignment file $aln_file conteain more than one sequence with the same name... Pleas fix and rerun.\n";
		}
		if ($aln_file=~/([^\/]+)$/) {$gene_name=$1;}
		my $seq_name="";
		my $seq="";
		while (my $line=<ALN>)
		{
			chomp ($line);
			if ($line=~/^>(.*)/)
			{
				# Process prev seq
				my $seq_length=length($seq);
				($Concate_Aln_Hash_ref,$Concate_Aln_Index_Hash_ref)=Add_seq_to_aln($seq_name,$seq,$gene_name,$total_aln_length,$Concate_Aln_Hash_ref,$Concate_Aln_Index_Hash_ref,$AlnLengths_ref) if ($seq_name ne "");
				
                # Start new seq
				$seq_name=$1; # new seq name
				$seq="";      # init
			}
			else
			{
				$seq=$seq.$line;
			}
		}
		# update last seq
		my $seq_length=length($seq);
		($Concate_Aln_Hash_ref,$Concate_Aln_Index_Hash_ref)=Add_seq_to_aln($seq_name,$seq,$gene_name,$total_aln_length,$Concate_Aln_Hash_ref,$Concate_Aln_Index_Hash_ref,$AlnLengths_ref);
		$total_aln_length=$total_aln_length+$seq_length if ($total_aln_length+$seq_length>$total_aln_length);
		# fill in missing 
		# print "fill_in_gaps($gene_name,$total_aln_length,$Concate_Aln_Hash_ref,$Concate_Aln_Index_Hash_ref);\n"; # DEBUG
		($Concate_Aln_Hash_ref,$Concate_Aln_Index_Hash_ref)=fill_in_gaps($gene_name,$total_aln_length,$Concate_Aln_Hash_ref,$Concate_Aln_Index_Hash_ref);

		# update the alignments lengths array
		push(@AlnLengths,$seq_length);
	}
}
close (LIST);
my %Info_Per_Pos=();
if ($PrintPerPos_Report eq "YES")
{
	open (LIST_OF_POS_DETAILS,$List_of_base_Aln_Info) || die "Can't open '$List_of_base_Aln_Info' $!";
	{
		while (my $Info_File=<LIST_OF_POS_DETAILS>)
		{
			print "## READING $Info_File";
			chomp ($Info_File);
			open (INFO_FILE,$Info_File) || die "Can't open '$Info_File' $!";
			my $line=<INFO_FILE>;
			if ($Info_File=~/(ENSG[0-9]+)/)
			{
				my $Gene_Name=$1;
				while (my $line=<INFO_FILE>)
				{
					chomp ($line);
					my ($Pos_On_Base_Aln,$Species_code,$Indel_Start_Relative_To_Seq,$Indel_End_Relative_To_Seq,$Alignments_Agree)=split(/\t/,$line);
					if (exists $Info_Per_Pos{$Gene_Name}{$Pos_On_Base_Aln})
					{
						$Info_Per_Pos{$Gene_Name}{$Pos_On_Base_Aln}=$Info_Per_Pos{$Gene_Name}{$Pos_On_Base_Aln}.";$Species_code,$Indel_Start_Relative_To_Seq,$Indel_End_Relative_To_Seq,$Alignments_Agree";
					}
					else
					{
						$Info_Per_Pos{$Gene_Name}{$Pos_On_Base_Aln}="$Species_code,$Indel_Start_Relative_To_Seq,$Indel_End_Relative_To_Seq,$Alignments_Agree";
					}
				}
			}
			close (INFO_FILE);

		}
	}
	close (LIST_OF_POS_DETAILS);
}
open (OUT_ALN,">$OutAln") || die "Can't open Concatenated Alignment for writing: '$OutAln' $!";
foreach my $species (keys %{$Concate_Aln_Hash_ref})
{
	print OUT_ALN ">$species\n$Concate_Aln_Hash_ref->{$species}\n";
}
close OUT_ALN;
open (OUT_REPORT,">$Concate_report") || die "Can't open the report '$Concate_report' for writing: $!";
print OUT_REPORT "species\tstart_on_concate_aln\tend_on_concate_aln\tgene_name\n";
open (OUT_PER_POS_REPORT,">$Concate_report_Per_Pos") || die "Can't open the report '$Concate_report_Per_Pos' for writing: $!" if ($PrintPerPos_Report eq "YES");
print OUT_PER_POS_REPORT "Pos_on_Concat_Aln\tSpecies\tIndel_Start_Relative_To_Seq\tIndel_End_Relative_To_Seq\tgene_name\tPos_On_Base_Aln\tAlignments_Agree\n" if ($PrintPerPos_Report eq "YES");
#Pos_on_Concat_Aln       Species Start_on_SeqBase_Aln    End_on_Base     Gene    Pos_On_Base_Aln Reliaibility
#Pan     0       16      ENSG00000000005_TNMD_sh_AA_MSA.MAFFT.SIC_CODED_WITH_REAL_NAMES.mfa
#species,
my $species_processed=0; # TO CONSIDER CHANGE
foreach my $species (keys %{$Concate_Aln_Index_Hash_ref})
{
	$species_processed++;
	foreach my $gene_index (@{$Concate_Aln_Index_Hash_ref->{$species}})
	{
		my ($start_on_concate_aln,$end_on_concate_aln,$gene_name)=split(",",$gene_index);
		print OUT_REPORT "$species\t$start_on_concate_aln\t$end_on_concate_aln\t$gene_name\n";
		if (($PrintPerPos_Report eq "YES") and ($species_processed<2))
		{
			if ($gene_name=~/(ENSG[0-9]+)/)
			{
				$gene_name=$1;
			}
			my $length_of_base_aln=$end_on_concate_aln-$start_on_concate_aln; 
			my $Pos_on_Concat_Aln=$start_on_concate_aln+1; # start_on_concate_aln starts from 0
			for (my $Pos_On_Base_Aln=1;$Pos_On_Base_Aln<=$length_of_base_aln;$Pos_On_Base_Aln++)
			{
				if (defined $Info_Per_Pos{$gene_name}{$Pos_On_Base_Aln})
				{
					my @Indel_Info=split(";",$Info_Per_Pos{$gene_name}{$Pos_On_Base_Aln});
					foreach my $Indel_Precense (@Indel_Info)
					{
						my ($Species_code,$Indel_Start_Relative_To_Seq,$Indel_End_Relative_To_Seq,$Alignments_Agree)=split(",",$Indel_Precense);
						my $Species="";
						if ($Seq_codes_to_Names ne "NA") {$Species=$Seq_codes_to_Names{$Species_code};}
						else {$Species=$Species_code;}
						
						print OUT_PER_POS_REPORT "$Pos_on_Concat_Aln\t$Species\t$Indel_Start_Relative_To_Seq\t$Indel_End_Relative_To_Seq\t$gene_name\t$Pos_On_Base_Aln\t$Alignments_Agree\n";
					}
					$Pos_on_Concat_Aln++;
				}
				else
				{
					print "[WARNNING] Info_Per_Pos{$gene_name}{$Pos_On_Base_Aln} - UNDEFINED\n";
				}
			}
		}
	}
}
close (OUT_REPORT);

sub Add_seq_to_aln
{
	my $seq_name=shift;          # seq name to add
	my $seq=shift;               # seq to add
	my $gene_name=shift;         # The aln the seq came from
	my $concate_aln_length=shift;# The total length of the concatenated aln
	my $Aln_Hash_ref=shift;      # ref to aln hash (key: seq_name (expected to be species name), value: concatenated seq)
	my $Aln_Data_ref=shift;      # ref to hash (key: species; value: array of start,end,Aln_name)
	my $Alignments_Length_ref=shift; # ref to ana array in which each elemant holt the Ith alignment length
	my $seq_length=length($seq);
	if (exists $Aln_Hash_ref->{$seq_name}) # Already was on the Aln
	{
		if (length ($Aln_Hash_ref->{$seq_name})==$concate_aln_length) # the species appeared on all prev instances
		{
			$Aln_Hash_ref->{$seq_name}=$Aln_Hash_ref->{$seq_name}.$seq;
			my $start=$concate_aln_length;
			my $end=$concate_aln_length+$seq_length;
			push (@{$Aln_Data_ref->{$seq_name}},"$start,$end,$gene_name");
			print "ADD [Already was on the Aln]: $seq_name\t$start\t$end\t$gene_name\n" if ($Debug eq "YES");
		}
		else # the species was missing only for prev gene
		{
			my $species_aln_length=length($Aln_Hash_ref->{$seq_name});
			my $fill_in_size=$concate_aln_length-$species_aln_length;
			my $fill_in_str="?" x $fill_in_size;
			$Aln_Hash_ref->{$seq_name}=$Aln_Hash_ref->{$seq_name}.$fill_in_str.$seq;
			print "[WARNNING Add_seq_to_aln] speicies: '$seq_name' is missing for gene: '$gene_name'; filling in with \"?\" from $species_aln_length to $concate_aln_length\n" if ($Debug eq "YES");
			my $TheNumOfPrevAln=scalar(@{$Alignments_Length_ref});
			if (exists $Aln_Data_ref->{$seq_name}) {push (@{$Aln_Data_ref->{$seq_name}},"$species_aln_length,$concate_aln_length,? [$TheNumOfPrevAln]");}
			# NOT RELEVANT AS WE ARE ALREADY IN UPPER IF THAT CHECK THE: if (exists $Aln_Hash_ref->{$seq_name}
#			else {
#				# $Aln_Data_ref->{$seq_name}=["$species_aln_length,$concate_aln_length,?"]; # was before added the for above
#				# on the data hash we want to account for all options it was missing before and create line for each missing
#				my $cummSum=0;
#				for (my $i=1;$i<scalar(@{$Alignments_Length_ref});$i++)
#				{
#					$cummSum=$cummSum+$Alignments_Length_ref->[$i];
#					print "QA: in foreach2: ",$cummSum-$Alignments_Length_ref->[$i].",".$cummSum.",?\n";
#					if ($cummSum<$end) 
#					{
#						my $LocalStart=$cummSum-$Alignments_Length_ref->[$i];
#						push (@{$Aln_Data_ref->{$seq_name}},"$LocalStart,$cummSum,? [$i]");
#						print "ADD [the species was missing only for prev gene]:$seq_name\t".$LocalStart.",".$cummSum.",? [$i]\n";
#					}
#				}
#			}
#			
#
#
##	for (my $i=1;$i<scalar(@{$Alignments_Length_ref});$i++)
##				{
##					print "QA: in foreach1: ".$Alignments_Length_ref->[$i-1].",".$Alignments_Length_ref->[$i].",?\n";
##					if ($Alignments_Length_ref->[$i]<$species_aln_length)  ### HERE 
##					{
##						$Aln_Data_ref->{$seq_name}=[$Alignments_Length_ref->[$i-1].",".$Alignments_Length_ref->[$i].",?"];
##					}
##				}
##                #$Aln_Data_ref->{$seq_name}=["$species_aln_length,$concate_aln_length,?"]; # WAS BEFORE THE FOREACH
##			}
#
#
#
##			push (@{$Aln_Data_ref->{$seq_name}},$concate_aln_length.",".$concate_aln_length+$seq_length.",".$gene_name);                                   # was before added the for above
##			print "ADD [the species was missing only for prev gene]: $seq_name\t$concate_aln_length\t".$concate_aln_length+$seq_length."\t$gene_name\n";   # was before added the for above 
		}
	}
	
	else # the species was not seen before
	{
		if ($concate_aln_length==0) # first gene
		{
			$Aln_Hash_ref->{$seq_name}=$seq;
			$Aln_Data_ref->{$seq_name}=["0,$seq_length,$gene_name"];
			print "ADD [first gene]:\t$seq_name\t0\t$seq_length\t$gene_name\n" if ($Debug eq "YES");
		}
		else # missing in all gene before
		{
			$Aln_Hash_ref->{$seq_name}=$seq;
			my $fill_in="?" x $concate_aln_length;
			$Aln_Hash_ref->{$seq_name}=$fill_in.$seq;
			print "[WARNNING Add_seq_to_aln] Organism: '$seq_name' First appeared on on Gene: $gene_name, Filled $concate_aln_length of \"?\"\n" if ($Debug eq "YES");

#			push (@{$Aln_Data_ref->{$seq_name}},"0,$concate_aln_length,?");                              # was used before added the for below
#			print "ADD [missing in all gene before]:$seq_name\t0\t".$concate_aln_length."\t?\n";         # was used before added the for below

#			push (@{$Aln_Data_ref->{$seq_name}},"$concate_aln_length,".$concate_aln_length+$seq_length.",$gene_name"); 
			my $end=$concate_aln_length+$seq_length;

			# on the data hash we want to account for all options it was missing before and create line for each missing
			my $cummSum=0;
			for (my $i=1;$i<scalar(@{$Alignments_Length_ref});$i++)
			{
				$cummSum=$cummSum+$Alignments_Length_ref->[$i];
				print "QA: in foreach2: ",$cummSum-$Alignments_Length_ref->[$i].",".$cummSum.",?\n" if ($Debug eq "YES");
				if ($cummSum<$end)  ### HERE 
				{
					my $LocalStart=$cummSum-$Alignments_Length_ref->[$i];
					push (@{$Aln_Data_ref->{$seq_name}},"$LocalStart,$cummSum,? [$i]");
					print "ADD [missing in all gene before]:$seq_name\t".$LocalStart.",".$cummSum.",? [$i]\n" if ($Debug eq "YES");
				}
			}
#			push (@{$Aln_Data_ref->{$seq_name}},"$concate_aln_length,$end,$gene_name");                    # was before added the for
#			print "ADD [missing in all gene before]:$seq_name\t$concate_aln_length\t$end\t$gene_name\n";   # was before added the for
			
			# add this seq
			push (@{$Aln_Data_ref->{$seq_name}},"$cummSum,$end,$gene_name");                   
			print "ADD [missing in all gene before]:$seq_name\t$cummSum\t$end\t$gene_name\n" if ($Debug eq "YES");  
			
#			print "AFTER ADD: $end\n";
		}
	}
	return ($Aln_Hash_ref,$Aln_Data_ref);
}
sub fill_in_gaps
{
	my $gene_name=shift;          # The aln the seq came from
	my $concate_aln_length=shift; # The total length of the concatenated aln
	my $Aln_Hash_ref=shift;       # ref to aln hash (key: seq_name (expected to be species name), value: concatenated seq)
	my $Aln_Data_ref=shift;       # ref to hash (key: species; value: array of start,end,Aln_name)

	foreach my $species (keys %{$Aln_Hash_ref})
	{
		# print "SPECIES:$species\n"; # DEBUG
		my $aln_length=length($Aln_Hash_ref->{$species});
		if ($concate_aln_length!=$aln_length)
		{
			my $fill_in_size=$concate_aln_length-$aln_length;
			my $fill_in_str="?" x $fill_in_size;
			$Aln_Hash_ref->{$species}=$Aln_Hash_ref->{$species}.$fill_in_str;
			print "[WARNNING fill_in_gaps] speicies: '$species' is missing for gene: '$gene_name'; filling in with \"?\" from $aln_length to $concate_aln_length\n" if ($Debug eq "YES");
			if (exists $Aln_Data_ref->{$species}) {push (@{$Aln_Data_ref->{$species}},"$aln_length,$concate_aln_length,? [for $gene_name]");}
			else {$Aln_Data_ref->{$species}=["$aln_length,$concate_aln_length,? [for $gene_name]"];}
		}
	}
	return ($Aln_Hash_ref,$Aln_Data_ref);
}
sub validate_seq_names_unique
{
	my $fasta_file=shift;
	my %unique_hash_names=();
	open (FASTA,$fasta_file)||die "Can't open FASTA file: '$fasta_file' $!";
	while (my $line=<FASTA>)
	{
		chomp ($line);
		if ($line=~/^>(.*)/)
		{
			my $seq_name=$1;
			if (exists $unique_hash_names{$seq_name})
			{
				$unique_hash_names{$seq_name}++;
			}
			else
			{
				$unique_hash_names{$seq_name}=1;
			}
		}
	}
	close (FASTA);
	my @duplicates=();
	foreach my $seq_name (keys %unique_hash_names)
	{
		if ($unique_hash_names{$seq_name}>1)
		{
			push (@duplicates,$seq_name);
		}
	}
	if (scalar(@duplicates)>0)
	{
		return ("ERROR",\@duplicates);
	}
	else
	{
		return ("OK");
	}
}
