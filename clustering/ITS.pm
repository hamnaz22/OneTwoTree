package ITS;
use File::Copy;
use Bio::SeqIO;
#use lib '/groups/itay_mayrose/annarice/ploidb/Miscellaneous';
use OneSeq;
use List::Util qw( min max );

sub pickOneITSTypePerSpeciesFasta{
	if(@_ < 6){
		die "usage:
		(1) input file - FASTA
		(2) genBank index
		(3) output directory
		(4) output file
		(5) scripts directory
		(6) log file
		";
	}


	my ($inputFile, $gbIndex, $outDir, $outputFile, $scriptsDir, $logFile) = @_;
	open(LOG,">$logFile") or die "can't open file $logFile";
	my ($fastaSeqObj, $header, $chosenSeq, $seq1Obj, $seq2Obj, $speciesDir, $gi, $append);
	my $c=0; #DEBUG

    print ("Calling1 perl $scriptsDir/SplitFastaSeqsBySpecies.pl $inputFile $outDir/species_all\n");
	system "perl $scriptsDir/SplitFastaSeqsBySpecies.pl $inputFile $outDir/species_all";

	my @chosenSequences = ();

	opendir( SP, "$outDir/species_all") or die "can't open $outDir/species_all";
	my @speciesFiles = grep { $_ ne '.' && $_ ne '..' } readdir(SP);
	my $out = Bio::SeqIO->new("-file" => ">$outputFile", "-format" => "Fasta");

    foreach my $species (@speciesFiles){
        print LOG "Processing $species\n";

        my $voucherUsed = 0;
        if($species =~ m/^(.+)\.fasta/){
            $speciesDir = $1;
        }
        else{
            next;
        }

        $speciesDir = "$outDir/species_all/$speciesDir";
        mkdir ($speciesDir) unless (-d $speciesDir);
        my ($ITS1count, $ITS2count, $combinedCount) = splitITS("$outDir/species_all/$species", "$speciesDir/ITS1_only", "$speciesDir/ITS2_only", "$speciesDir/combined", $gbIndex);

        # pick a combined sequence if available:
        if($combinedCount > 0){
            $chosenSeq = pickFromFile($combinedCount, "$speciesDir/combined", $speciesDir);
            push(@chosenSequences,$chosenSeq);
            print LOG "species $species: combined\t".getGI($chosenSeq)." was chosen\n"; #DEBUG
        }

        # else, merge two separated sequences - ITS1 + ITS2:
        elsif($ITS1count > 0 || $ITS2count > 0){

            #if there is voucher info:
            my ($isVoucherITS1, $hashRef1, $isVoucherITS2, $hashRef2) = (0,0,0,0);
            my (@keys1, @keys2, @appended);

            #in case both ITS1 and ITS2 exist, check if voucher information is available:
            if($ITS1count > 0 && $ITS2count > 0){
                ($isVoucherITS1, $hashRef1) = voucher("$speciesDir/ITS1_only", $gbIndex);
                ($isVoucherITS2, $hashRef2) = voucher("$speciesDir/ITS2_only", $gbIndex);

                if($isVoucherITS1 && $isVoucherITS2){
                    @keys1 = keys(%{$hashRef1}); #get the vouchers appear in ITS1
                    @keys2 = keys(%{$hashRef2}); #get the vouchers appear in ITS2

                    foreach my $voucherID (@keys1){ # append two sequences from the same voucherID
                        if(exists $hashRef2->{$voucherID}){
                            # pick ITS1 and ITS2 sequences with the same voucher ID, and then put them together
                            # (merge or append - according to whether they come from MSA or not)
                            $seq2Obj = pickFromArr($hashRef2->{$voucherID}, "abc2", $speciesDir);
                            $seq1Obj = pickFromArr($hashRef1->{$voucherID}, "abc1", $speciesDir);

                            push(@chosenSequences,$seq1Obj);
                            push(@chosenSequences,$seq2Obj);
                            $voucherUsed = 1;
                            # TODO MOSHE: Do we want just to take the first one with voucher?
                            last;
                        }
                    }
                }
            }

            if(!$voucherUsed){ #no use of voucher ID
                if($ITS1count > 0) {
                    $seq1Obj = pickFromFile($ITS1count, "$speciesDir/ITS1_only", $speciesDir);
                    push(@chosenSequences,$seq1Obj);
                }
                if($ITS2count > 0) {
                    $seq2Obj = pickFromFile($ITS2count, "$speciesDir/ITS2_only", $speciesDir);
                    push(@chosenSequences,$seq2Obj);
                }
            }
        }

        $c++; #DEBUG
    }

    my $numOfChosen = @chosenSequences;
    print ("A total of $numOfChosen seqs were chosen\n");
    foreach (@chosenSequences) {
        my $seqToWrite = $_;
        print ("Writing seq ".getGI($seqToWrite)."\n");
        $out->write_seq($seqToWrite);
    }

    print LOG "Number of species = $c\n";
    close(LOG);
}

sub pickOneSeqPerSpeciesMSA{
	if(@_ < 7){
		die "usage:
		(1) input file - MSA
		(2) genBank index
		(3) output directory
		(4) output file
		(5) scripts directory
		(6) log file
		(7) append
		";
	}

	my ($inputFile, $gbIndex, $outDir, $outputFile, $scriptsDir, $logFile,$append) = @_;
	open(LOG,">$logFile") or die "can't open file $logFile";
	my ($fastaSeqObj, $header, $chosenSeq, $seq1Obj, $seq2Obj, $speciesDir, $gi);
	my $c=0; #DEBUG

    print ("Calling2 perl $scriptsDir/SplitFastaSeqsBySpecies.pl $inputFile $outDir/species\n");
	system "perl $scriptsDir/SplitFastaSeqsBySpecies.pl $inputFile $outDir/species";
	
	my @chosenSequences = ();
	
	opendir( SP, "$outDir/species") or die "can't open $outDir/species";
	my @speciesFiles = grep { $_ ne '.' && $_ ne '..' } readdir(SP);
	my $out = Bio::SeqIO->new("-file" => ">$outputFile", "-format" => "Fasta");

    foreach my $species (@speciesFiles){
        print LOG "Processing $species\n";

        if($species =~ m/^(.+)\.fasta/){
            $speciesDir = $1;
        }
        else{
            next;
        }

        $speciesDir = "$outDir/species/$speciesDir";
        mkdir ($speciesDir) unless (-d $speciesDir);
        my ($ITS1count, $ITS2count, $combinedCount) = splitITS("$outDir/species/$species", "$speciesDir/ITS1_only", "$speciesDir/ITS2_only", "$speciesDir/combined", $gbIndex);

        if ($ITS1count > 1 || $ITS2count > 1 || $combinedCount > 1) {
            die "More than a single ITS per type was found in $speciesDir $ITS1count $ITS2count $combinedCount";
        }
        if (($ITS1count > 1 || $ITS2count > 1) && ($combinedCount > 0)) {
            die "Combined ITS was found as well as ITS1/ITS2 $ITS1count $ITS2count $combinedCount";
        }

        if($combinedCount > 0){
            # pick a combined sequence if available:
            $chosenSeq = pickFromFile($combinedCount, "$speciesDir/combined", $speciesDir);
            print LOG "species $species: combined\t".getGI($chosenSeq)." was chosen\n"; #DEBUG
        } elsif($ITS1count > 0 && $ITS2count > 0){
            # else, merge two separated sequences - ITS1 + ITS2:
            $seqObj1 = getTheFirstSequence("$speciesDir/ITS1_only");
            $seqObj2 = getTheFirstSequence("$speciesDir/ITS2_only");
            print LOG "Merging/appending ITS1 ".getGI($seqObj1)." and ITS2 ".getGI($seqObj2)."\n"; #DEBUG
            $chosenSeq = getChosenSequence(1, 1, $seqObj1, $seqObj2, $append);
            print LOG "species $species: ITS1 and ITS2 were merged \t".getGI($chosenSeq)."\n"; #DEBUG
        } elsif($ITS1count > 0) {
            print LOG "species $species: using ITS1";
            $chosenSeq = getTheFirstSequence("$speciesDir/ITS1_only\n");
        } elsif($ITS2count > 0) {
            print LOG "species $species: using ITS2";
            $chosenSeq = getTheFirstSequence("$speciesDir/ITS2_only\n");
        } else {
            print STDERR "ERROR - species $species: has no ITS sequence \n";
        }

        $out->write_seq($chosenSeq);
        $c++; #DEBUG
    }
    print LOG "Number of species = $c\n";
    close(LOG);
}

sub getGI{
	my ($seqObj) = @_;
	#my ($GI, $t) = ($seqObj->id() =~ m/gi\|(\d+)\|taxonid\|(\d+)\|/); 
	my ($GI, $t) = ($seqObj->id() =~ m/gi\|([A-Z{2}\d]+\.\d+)\|taxonid\|(\d+)\|/);
	return "$GI";
}


=begin
sub fastaFileToSeqObjArray{
	my ($file) = @_;
	my $fastaIn = Bio::SeqIO->new("-format" => 'fasta', "-file"  => $file);
	my ($seqObj, @seqObjArray);
	while( $seqObj = $fastaIn->next_seq() ) {
		push(@seqObjArray,$seqObj);
	}
	return @seqObjArray;
}
=cut

sub fastaFileToNumberOfSeqs{
	my ($file) = @_;
	my $fastaIn = Bio::SeqIO->new("-format" => 'fasta', "-file"  => $file);
	my $count = 0;
	while( $seqObj = $fastaIn->next_seq() ) {
        $count++;
	}
	return $count;
}


sub createGIArrayFromSeqObj{
	my (@seqObjArr) = @_;
	my @giArr;
	foreach my $seqObj (@seqObjArr){	
		#if ($seqObj->id() =~ m/gi\|(\d+)\|/){
		if ($seqObj->id() =~ m/gi\|([A-Z{2}\d]+\.\d+)\|/){
			push(@giArr, $1);
		}
	}
	return @giArr;
}

sub pickSeq{ ### call Anna's script
	my ($inFileFasta, $outDir) = @_;
	open(IN, "<$inFileFasta") or die "can't open file $inFileFasta";
	my @fastaLines = <IN>;
	if(@fastaLines > 0){	
		#generate all-v-all blast:
		my ($blast) = ($inFileFasta =~ m/.*\/([^\/]+)$/);
		$blast = "$outDir/$blast"."_all_v_all.blastn";
		system "formatdb -i $inFileFasta -pF";    #-pF -- nucleotides
		system "blastall -p blastn -d $inFileFasta -i $inFileFasta -v 100000 -b 100000 -e 1e-5 -m 8 > $blast";
		
		#get the hash of the blast:
		my ($de_combinationScores,$min)=OneSeq::BuildBlastHash($blast,$inFileFasta);
		
		#get the selected gi:
		my $giArrRef = OneSeq::CreateGIArray([@fastaLines]);
		print("!!!!!!!!!!!!!!!!!!!ReturnMaxBitScore ReturnMaxBitScore ReturnMaxBitScore params:\n");
		print ("CreateGIArray return valuse:\n");
		foreach(@giArrRef)
		{
			print "$_\r\n";
		}
		
		print ("BuildBlastHash return valuse:  $de_combinationScores \n");
		print ("min return valuse:  $min \n");
		my $selectedGI = OneSeq::ReturnMaxBitScore($giArrRef , $de_combinationScores, $min, $inFileFasta); #seqfile and cluster - undef (temporary)
		#get the sequence object of this gi:
		$in = Bio::SeqIO->new(-format=>'fasta', -file=>$inFileFasta);
		my $seqObj;
		while (defined($seqObj = $in->next_seq())){
			if (($seqObj->id()) =~ m/$selectedGI/){
				return $seqObj;
			}
		}	
	}
	else{
		die "pickSeq - empty array!\n";
	}
	
	close(IN);
}

sub pickFromFile{
	my ($count, $file, $outDir)=@_;
	my $seqObj;
	if($count > 1){
		$seqObj = pickSeq($file, $outDir);
	}
	elsif($count == 1){
		$seqObj = getTheFirstSequence($file); 
	}
	return $seqObj;
}

#in case of array, write it to file and send to pickFromFile
sub pickFromArr{
	my ($seqObjArrRef, $outFile, $outDir)=@_;
	@seqObjArr = @{$seqObjArrRef};
	if(@seqObjArr == 1){
		return $seqObjArr[0];
	}
	else{
		my $out = Bio::SeqIO->new("-file"=>">$outFile","-format"=>"Fasta");
		foreach my $seqObj (@seqObjArr){
			$out->write_seq($seqObj);
			#$seqObj->write_seq($out);
		}
		return pickFromFile(scalar(@seqObjArr), $outFile, $outDir);
	}
}

sub getChosenSequence{
	
	my ($ITS1count, $ITS2count, $seq1Obj, $seq2Obj, $append) = @_;
	my $chosenSeq;
	if($ITS1count > 0 && $ITS2count > 0){
		if($append){
			$chosenSeq = appendSequences($seq1Obj, $seq2Obj);				
		}
		else{
			$chosenSeq = mergeSequences($seq1Obj, $seq2Obj);				
		}
	}
	elsif($ITS1count > 0){
		$chosenSeq = $seq1Obj;
	}
	elsif($ITS2count > 0){
		$chosenSeq =  $seq2Obj;
	}
	return $chosenSeq;
}

sub voucher{
	
	my ($inFasta, $gbIndex) = @_;
	my %hash = (); #hash to contain voucher-->sequence object(s) of this voucher
	my (@sourceValues, @features, $gbSeqObj,$GI, $voucherID);	
	
	my $in = Bio::SeqIO->new("-file"=>"<$inFasta","-format"=>"Fasta") or die "can't open $inFasta";	#open all-seqs fasta file

	while (defined($fastaSeqObj = $in->next_seq())) { #iterate over all the sequences
		$header = getHeader($fastaSeqObj); #get the header of the fasta entry
		#($GI) = ($header =~ m/gi\|(\d+)\|/) ;
		($GI) = ($header =~ m/gi\|([A-Z{2}\d]+\.\d+)\|/) ;

		$gbSeqObj = $gbIndex->fetch($GI); # get the entry from GenBank

		@features = $gbSeqObj->get_SeqFeatures; #get the sequence's features (cds, source, gene....)

		foreach my $feature (@features){
			if ($feature->primary_tag eq "source"){
				@sourceValues = $feature->get_all_tags;
				foreach my $tag (@sourceValues){
					if ($tag eq "isolate"){
						$voucherID = join( ' ', $feature->get_tag_values($tag));
						goto label;
					}
				}	
				
				foreach my $tag (@sourceValues){
					if ($tag eq "specimen_voucher"){
						$voucherID = join( ' ', $feature->get_tag_values($tag));
					}
				}
				
				label:
				if(defined $voucherID){
#					my $id = $fastaSeqObj->id(); #DEBUG
#					my $desc = $fastaSeqObj->desc();#DEBUG
#					my $seq = $fastaSeqObj->seq();#DEBUG
					push(@{$hash{$voucherID}}, $fastaSeqObj);
				}
			}
		}
	}
	
	return (defined $voucherID ? 1 : 0), {%hash};

}

sub getTheFirstSequence{
	if (@_ < 1){
		die "usage: (1) input file - MSA in fasta format";
	}
	my ($inFile) = @_;
	#print "Getting first seq from $inFile\n";
	my $in = Bio::SeqIO->new("-file"=>"<$inFile","-format"=>"Fasta");	#open all-seqs fasta file
	my $firstSeq = $in->next_seq();
	#print "Returning $firstSeq as first seq in $inFile\n";
	return $firstSeq;
}

#append two given sequences (not aligned to each other)
sub appendSequences{
	my ($seq1Obj, $seq2Obj) = @_;
	my $appendedSeq = $seq1Obj->seq() . $seq2Obj->seq();
	
	my $appendedSeqObj = createMergedSeqObj($seq1Obj, $seq2Obj, $appendedSeq);
		
	return $appendedSeqObj;

}

#create sequence object for the merged sequence created from seq1Obj and seq2Obj
sub createMergedSeqObj{
	my ($seq1Obj, $seq2Obj, $mergedSeq) = @_;
	my $header1 = getHeader($seq1Obj);
	my $header2 = getHeader($seq2Obj);
	my ($gi1) = ($header1 =~ m/gi\|([^\|]+)\|/);
	my ($gi2) = ($header2 =~ m/gi\|([^\|]+)\|/);
	my $newHeader = "$header1|ITS-merge|ITS1: $gi1 ITS2: $gi2";
	my ($id, $desc) = ($newHeader =~ /([^ ]*) (.*)/); #split the header into id and description (by the bioperl's genbank format): id = the header until the first space, desc = the rest of the header
	my $mergedSeqObj = Bio::Seq->new(-id => $id, -desc  => $desc, -seq => $mergedSeq);
	return $mergedSeqObj;
}

# merge two given aligned sequences (aligned to each other)
sub mergeSequences{

	my ($seq1Obj, $seq2Obj) = @_;
	my (@mergedSeq);
	my @seq1 = split(//, $seq1Obj->seq());
	my @seq2 = split(//, $seq2Obj->seq());
		
	for(my $i=0; $i< min(scalar(@seq1), scalar(@seq2)); $i++){
		if($seq1[$i] eq "-" && $seq2[$i] ne "-"){
			$mergedSeq[$i] = $seq2[$i];
		}
		else{
			$mergedSeq[$i] = $seq1[$i];
		}
	}
	
	my $mergedSeqObj = createMergedSeqObj($seq1Obj, $seq2Obj, join("",@mergedSeq));		
	return $mergedSeqObj;

}

# create MSA of the ITS sequences
sub ITS_MSA{
	if (@_ < 8){
		die "usage:
		(1) output directory
		(2) number of ITS1 sequences
		(3) number of ITS2 sequences
		(4) number of combined sequences
		(5) fasta file of ITS1 sequences
		(6) fasta file of ITS2 sequences
		(7) fasta file of combined sequences
		(8) external scripts dir
		";
	}
	my ($outDir, $ITS1count, $ITS2count, $combinedCount, $ITS1Fasta, $ITS2Fasta, $combinedFasta,$scriptsDir) = @_;
	my ($outputFileName, $combined);
	
	#combined sequences exist
    print ("About to adjust direction for its sequences\n");
	if ($combinedCount > 0){
	    # Adjusting the direction of the relevant fasta files
        if($ITS1count > 0 || $ITS2count > 0){
            system "cat $ITS1Fasta > $outDir/SEP_ITS1+ITS2.fasta";
            system "cat $ITS2Fasta >> $outDir/SEP_ITS1+ITS2.fasta";

            print ("Calling $scriptsDir/adjustdir.sh $combinedFasta,$outDir/SEP_ITS1+ITS2.fasta\n");	# Michal: Changed path $scriptsDir/../ploidb/shell/adjustdir.sh
            system "$scriptsDir/adjustdir.sh $combinedFasta,$outDir/SEP_ITS1+ITS2.fasta";
        } else {
            print ("Calling $scriptsDir/adjustdir.sh $combinedFasta\n");
            system "$scriptsDir/adjustdir.sh $combinedFasta";
        }
		
		#create MSA of the combined sequences (contain both ITS1 and ITS2)
		if ($combinedCount > 1){
			$combined = "$outDir/combined.msa";
			#system "mafft --retree 2 --maxiterate 0 --bl 62 --op 1.530000 --ep 0.000000 $combinedFasta > $combined" ;
			system "mafft --auto --ep 0.000000 $combinedFasta > $combined" ;
		}
		else {
			if ($combinedCount eq 1){
				$combined = "$outDir/combined.msa";
				open($IN, "$combinedFasta" ) or die "ERROR:cannot open file $!";
				open($OUT, ">$combined") or die "ERROR:cannot create file $!"; 
				my @totalcount = <$IN>;
				foreach(@totalcount){
				print $OUT $_;
				}
			}
			else{
				$combined = "combined.fasta";
			}	
		}			
		$outputFileName = $combined;
			
		#if there are separate sequences, add them to the MSA with 'addfragments' option of MAFFT
		if($ITS1count > 0 || $ITS2count > 0){
			system "mafft --addfragments $outDir/SEP_ITS1+ITS2.fasta --multipair $outDir/combined.msa > $outDir/combined+sep.msa";
			#system "mafft --addfragments $outDir/SEP_ITS1+ITS2.fasta --multipair $outDir/$combined > $outDir/combined+sep.msa";
			
			editITSdisplay("$outDir/combined+sep.msa"); #make a similar file, with the kind of ITS in the start of each sequence
			
			$outputFileName = "$outDir/combined+sep.msa";
		}
	}
	
	#no combined sequences
	else{  
	
		#MSA for ITS1 sequences (if more than one exist)
		if ($ITS1count > 1){
		    print ("Calling $scriptsDir/adjustdir.sh $ITS1Fasta\n");
            system "$scriptsDir/adjustdir.sh $ITS1Fasta";

			#system "mafft --retree 2 --maxiterate 0 --bl 62 --op 1.530000 --ep 0.000000 $ITS1Fasta > $outDir/ITS1_only.msa" ; #-bl 62 , --op 1.530000 -> is default, --ep 0.000000 allows large gaps !!!
			system "mafft --auto 0 --ep 0.000000 $ITS1Fasta > $outDir/ITS1_only.msa" ;
			$outputFileName =  "$outDir/ITS1_only.msa";
		}
		elsif($ITS1count == 1){
			$outputFileName = $ITS1Fasta;
		}
		
		#MSA for ITS2 sequences (if more than one exist)
		if ($ITS2count > 1){
		    print ("Calling $scriptsDir/adjustdir.sh $ITS2Fasta\n");
            system "$scriptsDir/adjustdir.sh $ITS2Fasta";

			#system "mafft --retree 2 --maxiterate 0 --bl 62 --op 1.530000 --ep 0.000000 $ITS2Fasta > $outDir/ITS2_only.msa" ;
			system "mafft --auto --ep 0.000000 $ITS2Fasta > $outDir/ITS2_only.msa" ;
			if(defined $outputFileName){
                system "cat $outputFileName > $outDir/SEP_ITS1+ITS2.msa";
                system "cat $outDir/ITS2_only.msa >> $outDir/SEP_ITS1+ITS2.msa";
				$outputFileName = "$outDir/SEP_ITS1+ITS2.msa"
			}
			else{
				$outputFileName = "$outDir/ITS2_only.msa";
			}
		}
		elsif($ITS2count == 1){
			if(defined $outputFileName){
                system "cat $outputFileName > $outDir/SEP_ITS1+ITS2.msa";
                system "cat $ITS2Fasta >> $outDir/SEP_ITS1+ITS2.msa";
				$outputFileName = "$outDir/SEP_ITS1+ITS2.msa"
			}
			else{
				$outputFileName = $ITS2Fasta;
			}
		}

		#Align SEP_ITS1_ITS2: added bu michal to ensure its aligned file
		system "cat $outDir/SEP_ITS1+ITS2.msa >> $outDir/SEP_ITS1+ITS2_temp.msa";
		system "mafft --auto --quiet SEP_ITS1+ITS2_temp.msa > $outputFileName" ;
	}
	
	# return the path of the actual output file
	return $outputFileName;
}

# create MSA of the ITS sequences
sub ITS_CLUSTALO{
	if (@_ < 8){
		die "usage:
		(1) output directory
		(2) number of ITS1 sequences
		(3) number of ITS2 sequences
		(4) number of combined sequences
		(5) fasta file of ITS1 sequences
		(6) fasta file of ITS2 sequences
		(7) fasta file of combined sequences
		(8) external scripts dir
		";
	}
	my ($outDir, $ITS1count, $ITS2count, $combinedCount, $ITS1Fasta, $ITS2Fasta, $combinedFasta,$scriptsDir) = @_;
	my ($outputFileName, $combined);

	#combined sequences exist
    print ("About to adjust direction for its sequences\n");
	if ($combinedCount > 0){
	    # Adjusting the direction of the relevant fasta files
        if($ITS1count > 0 || $ITS2count > 0){
            system "cat $ITS1Fasta > $outDir/SEP_ITS1+ITS2.fasta";
            system "cat $ITS2Fasta >> $outDir/SEP_ITS1+ITS2.fasta";

            print ("Calling $scriptsDir/adjustdir.sh $combinedFasta,$outDir/SEP_ITS1+ITS2.fasta\n");	# Michal: Changed path $scriptsDir/../ploidb/shell/adjustdir.sh
            system "$scriptsDir/adjustdir.sh $combinedFasta,$outDir/SEP_ITS1+ITS2.fasta";
        } else {
            print ("Calling $scriptsDir/adjustdir.sh $combinedFasta\n");
            system "$scriptsDir/adjustdir.sh $combinedFasta";
        }

		#create MSA of the combined sequences (contain both ITS1 and ITS2)
		if ($combinedCount > 1){
			$combined = "$outDir/combined.msa";
			#system "mafft --retree 2 --maxiterate 0 --bl 62 --op 1.530000 --ep 0.000000 $combinedFasta > $combined" ;
			#system "mafft --auto --ep 0.000000 $combinedFasta > $combined" ;
			system("/share/apps/clustalo -i $combinedFasta -o $combined --outfmt=fasta");
		}
		else {
			if ($combinedCount eq 1){
				$combined = "$outDir/combined.msa";
				open($IN, "$combinedFasta" ) or die "ERROR:cannot open file $!";
				open($OUT, ">$combined") or die "ERROR:cannot create file $!";
				my @totalcount = <$IN>;
				foreach(@totalcount){
				print $OUT $_;
				}
			}
			else{
				$combined = "combined.fasta";
			}
		}
		$outputFileName = $combined;

		#if there are separate sequences, add them to the MSA with 'addfragments' option of MAFFT
		if($ITS1count > 0 || $ITS2count > 0){
			system "/share/apps/clustalo --p2 $outDir/SEP_ITS1+ITS2.fasta --p1 $outDir/combined.msa -o $outDir/combined+sep.msa --outfmt=fasta"; #This option allows also 1 seq msa's to be combined
			#system "mafft --addfragments $outDir/SEP_ITS1+ITS2.fasta --multipair $outDir/combined.msa > $outDir/combined+sep.msa";
			#system "mafft --addfragments $outDir/SEP_ITS1+ITS2.fasta --multipair $outDir/$combined > $outDir/combined+sep.msa";

			editITSdisplay("$outDir/combined+sep.msa"); #make a similar file, with the kind of ITS in the start of each sequence

			$outputFileName = "$outDir/combined+sep.msa";
		}
	}

	#no combined sequences
	else{

		#MSA for ITS1 sequences (if more than one exist)
		if ($ITS1count > 1){
		    print ("Calling $scriptsDir/adjustdir.sh $ITS1Fasta\n");
            system "$scriptsDir/adjustdir.sh $ITS1Fasta";

			#system "mafft --retree 2 --maxiterate 0 --bl 62 --op 1.530000 --ep 0.000000 $ITS1Fasta > $outDir/ITS1_only.msa" ; #-bl 62 , --op 1.530000 -> is default, --ep 0.000000 allows large gaps !!!
			#system "mafft --auto 0 --ep 0.000000 $ITS1Fasta > $outDir/ITS1_only.msa" ;
			system "/share/apps/clustalo -i $ITS1Fasta -o $outDir/ITS1_only.msa --outfmt=fasta";
			$outputFileName =  "$outDir/ITS1_only.msa";
		}
		elsif($ITS1count == 1){
			$outputFileName = $ITS1Fasta;
		}

		#MSA for ITS2 sequences (if more than one exist)
		if ($ITS2count > 1){
		    print ("Calling $scriptsDir/adjustdir.sh $ITS2Fasta\n");
            system "$scriptsDir/adjustdir.sh $ITS2Fasta";

			#system "mafft --retree 2 --maxiterate 0 --bl 62 --op 1.530000 --ep 0.000000 $ITS2Fasta > $outDir/ITS2_only.msa" ;
			#system "mafft --auto --ep 0.000000 $ITS2Fasta > $outDir/ITS2_only.msa" ;
			system "/share/apps/clustalo -i $ITS2Fasta -o $outDir/ITS2_only.msa --outfmt=fasta";
			if(defined $outputFileName){
                system "cat $outputFileName > $outDir/SEP_ITS1+ITS2.msa";
                system "cat $outDir/ITS2_only.msa >> $outDir/SEP_ITS1+ITS2.msa";
				$outputFileName = "$outDir/SEP_ITS1+ITS2.msa"
			}
			else{
				$outputFileName = "$outDir/ITS2_only.msa";
			}
		}
		elsif($ITS2count == 1){
			if(defined $outputFileName){
                system "cat $outputFileName > $outDir/SEP_ITS1+ITS2.msa";
                system "cat $ITS2Fasta >> $outDir/SEP_ITS1+ITS2.msa";
				$outputFileName = "$outDir/SEP_ITS1+ITS2.msa"
			}
			else{
				$outputFileName = $ITS2Fasta;
			}
		}
	}

	# return the path of the actual output file
	return $outputFileName;
}


sub editITSdisplay{
	my ($finalFileName) = @_;
	
	if($finalFileName ne "" && $finalFileName !~ m/,/){
		open(IN, "<$finalFileName") or die "can't open file $finalFileName";
		open(OUT, ">$finalFileName.edited") or die "can't open file $finalFileName.edited";
		while (my $line = <IN>){
			$line =~ s/^>(.+)(\[(1|2|1\+2)\])/>$2$1/;
			print OUT $line;
		}
		close(IN);
		close(OUT);
	}
}

# given a fasta file, separate the sequences to 3 groups: (1) contain both ITS1 and ITS2, (2) contain ITS1 without ITS2, (3) contain ITS2 without ITS1 
sub splitITS{

	if(@_ < 5){
		die "usage:
		(1) input fasta file
		(2) output file contains all the sequences contain ITS1 and not ITS2
		(3) output file contains all the sequences contain ITS2 and not ITS1
		(4) output file contains all the sequences contain both ITS1 and ITS2
		(5) genBank index
		";
	}
	
	my ($input, $its1Output, $its2Output, $combinedOutput, $gbIndex) = @_;
	my ($fastaSeqObj, $formattedSeqObj, $header, $desc);
	my $in = Bio::SeqIO->new("-file"=>"<$input","-format"=>"Fasta");	#open all-seqs fasta file
	my $its1Out = Bio::SeqIO->new("-file" => ">$its1Output", "-format" => "Fasta");
	my $its2Out = Bio::SeqIO->new("-file" => ">$its2Output", "-format" => "Fasta");
	my $combinedOut = Bio::SeqIO->new("-file" => ">$combinedOutput", "-format" => "Fasta");
	
	
	my ($ITS1count, $ITS2count, $combinedCount) = (0,0,0);
	
	while (defined($fastaSeqObj = $in->next_seq())) { #iterate over all the sequences

		$header = getHeader($fastaSeqObj);
		#my ($GI) = ($header =~ m/gi\|(\d+)\|/);
		my ($GI) = ($header =~ m/gi\|([A-Z{2}\d]+\.\d+)\|/);
		$desc = getFeaturesInfo($GI, $gbIndex);
		$desc = formatKeyWords($desc);

		#format the header to remove unnecessary ">" characters
		$header =~ s/<(.+)>/[$1]/g;
		$header =~ s/(>|<)//g;
		my ($a, $b) = ($header =~ m/^(.+\|description\|)(.+)$/);

		if($desc =~ m/ITS1/ && $desc =~ m/ITS2/){
			#$header = "$a"."[1+2]"."$b"; #mark the sequence as a combined one
			my ($id, $desc1) = ($header =~ /([^ ]*) (.*)/); #split the header into id and description (by the bioperl's genbank format): id = the header until the first space, desc = the rest of the header
			$formattedSeqObj = Bio::Seq->new(-id => $id, -desc  => $desc1, -seq => $fastaSeqObj->seq() );
			
			$combinedOut->write_seq($formattedSeqObj);
			$combinedCount++;
		}
		elsif($desc =~ m/ITS1/ && $desc !~ m/ITS2/){
			#$header = "$a"."[1]"."$b"; #mark the sequence as containing ITS1 only
			my ($id, $desc1) = ($header =~ /([^ ]*) (.*)/); #split the header into id and description (by the bioperl's genbank format): id = the header until the first space, desc = the rest of the header
			$formattedSeqObj = Bio::Seq->new(-id => $id, -desc  => $desc1, -seq => $fastaSeqObj->seq() );
			$its1Out->write_seq($formattedSeqObj);
			
			$ITS1count++;
		}
		elsif($desc !~ m/ITS1/ && $desc =~ m/ITS2/){
			#$header = "$a"."[2]"."$b"; #mark the sequence as containing ITS2 only
			my ($id, $desc1) = ($header =~ /([^ ]*) (.*)/); #split the header into id and description (by the bioperl's genbank format): id = the header until the first space, desc = the rest of the header
			$formattedSeqObj = Bio::Seq->new(-id => $id, -desc  => $desc1, -seq => $fastaSeqObj->seq() );
		
			$its2Out->write_seq($formattedSeqObj);
			$ITS2count++;
		}
	}

	return ($ITS1count, $ITS2count, $combinedCount);
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
# format the free-text to get a better collapse
sub formatKeyWords {
	my ($desc) = @_;
	$desc =~ s/internal trasncribed spacer 1/ITS1/ig;
	$desc =~ s/internal transcribed spacer 1/ITS1/ig;
	$desc =~ s/ITS 1/ITS1/ig;
	$desc =~ s/internal transcribed spacer 2/ITS2/ig;
	$desc =~ s/internal trasncribed spacer 2/ITS1/ig;
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

sub getHeader{
	my $fastaSeqObj = $_[0];
	my $id = $fastaSeqObj->id();
	my $desc = $fastaSeqObj->desc();
	my $header = $id." ".$desc;
	return $header;
}


# get the value of a given feature from the feature object
sub getFeatureValue {
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
		if ( $tag eq "mol_type" ) {
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

	if ( !defined $featureValue ) {
		$featureValue = $featureObject->primary_tag;
		return $featureValue;
	}

}




1;
