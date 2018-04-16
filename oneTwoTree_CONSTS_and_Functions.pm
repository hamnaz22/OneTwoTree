package oneTwoTree_CONSTS_and_Functions;

use strict;
use warnings;
use File::Slurp;
use File::Path;# qw(make_path);
use File::Basename;


use vars qw(@ISA @EXPORT);

require Exporter;

@ISA = qw(Exporter);

@EXPORT = qw(WriteToFile WriteToFileWithTimeStamp ReadFromFile CheckInputTypeValid ReadParam ReadToArrayParam CheckIfFileExists);



######### Consts ###########
use constant SERVER_NAME 				=> "oneTwoTree";
use constant RESULTS_DIR_ABSOLUTE_PATH 	=> "/bioseq/data/results/oneTwoTree";
use constant RESULTS_LINK 				=> "oneTwoTree_results";
use constant LOG_FILE 					=> "output.log";
use constant ERROR_STATUS_LOG_FILE 		=> "error.OU";
use constant OUTPUT_FILES_LIST 			=> "output_files.dat";
use constant LOG_DIR_ABSOLUTE_PATH 		=> "/bioseq/data/results/oneTwoTree/logs/";
use constant MAIN_PIPELINE_LOG 			=> "oneTwoTree.log";
#use constant PYTHON_MODULE_TO_LOAD 		=> 'python/python-3.3.0';
use constant PYTHON_MODULE_TO_LOAD 		=> 'python/anaconda3-5.0.0';
use constant PERL_MODULE_TO_LOAD 		=> 'perl/perl-5.16.3';
use constant QSUB_JOB_NUM_FILE			=> "qsub_job_num.dat";
use constant RESULTS_PAGE_URL			=> "/results.html";
use constant FILENAME_PARAMS_TXT		=> "params.txt";
use constant FILENAME_PARAMS_HTML		=> "params.html";
use constant FILENAME_LOCUS_TREENUM_TXT		=> "NumOfClusterTrees.txt";
use constant FILENAME_EMAIL				=> "email.txt";



######### Functions ############
sub WriteToFile
{
	my ($file, $message, $shouldOverwrite) = @_;
	
	# creating file dir if doesnt exist
	my ($fileName, $fileDir) = fileparse($file);
	#make_path($fileDir);
	mkpath($fileDir);
	
	$message =~ s/^\s+//;
	
	if (defined $shouldOverwrite && $shouldOverwrite){
		write_file($file, "$message\n");
	}
	else {
		append_file($file, "$message\n");
	}	
}


sub ReadParam
{
	my ($file, $paramName) = @_;
	my $NameValue = "NONE";
	#Check if file exists:
	if (-e $file)
	{
		open my $fh, '<', $file or die "Could not open '$file' $!\n";
		while (my $line = <$fh>) {
			chomp $line;
			if (index($line, $paramName) != -1)
			{
				my @fields = split /:/, $line;
				$NameValue = $fields[1];
				return $NameValue;	
			}
		}
	}
	$NameValue  = "NoData";
	return $NameValue;
}

sub ReadToArrayParam  
{
	my ($file, $paramName) = @_;
	my @values = ();
	
	#Check if file exists:
	if (-e $file)
	{
		open my $fh, '<', $file or die "Could not open '$file' $!\n";
		while (my $line = <$fh>) {
			chomp $line;
			@values = split /,/, $line;
			return \@values;
		}
	}
	return \@values;
}

sub CheckIfFileExists  
{
	my ($file) = @_;
	my $flag_exists = 'NO';
	
	#Check if file exists:
	if (-e $file)
	{
		$flag_exists = 'YES';
		return $flag_exists;
	}
	return $flag_exists;
}

sub ReadFromFile
{
	my ($file, $defaultValue) = @_;
	
	if (defined $defaultValue)
	{
		if (-e $file)
		{
			my $line = read_file($file);
			return $line;
		}
		else
		{
			# this is ugly. delete this after renaming all relevant files to dat.
			my ($name, $dir, $ext) = fileparse($file, ".dat");
		 	if ($ext eq ".dat")	{
		 		my $newFile = substr($file, 0, -4);
		 		
		 		if (-e $newFile) {
					my $line = read_file($newFile);
					return $line;
				}
		 	}
			
			return $defaultValue;
		}
	}
	else
	{
		my @lines;
		
		if (-e $file)
		{
			@lines = read_file($file);
		}
		else {
			# this is ugly. delete this after renaming all relevant files to dat.
			my ($name, $dir, $ext) = fileparse($file, ".dat");
		 	if ($ext eq ".dat")	{
		 		my $newFile = substr($file, 0, -4);
		 		
		 		if (-e $newFile) {
					@lines = read_file($newFile);
				}
		 	}
		}
		
		return @lines;
	}
}