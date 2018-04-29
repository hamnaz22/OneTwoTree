use strict;

if (@ARGV < 1){
	die "usage: input directory name";
}

my ( $inDirName ) = @ARGV;
my(@inFiles, $out, @lines, %filesHash);
%filesHash = ();
opendir( IN_DIR, "$inDirName" ) or die "can't open $inDirName";
@inFiles = grep { $_ ne '.' && $_ ne '..' } readdir(IN_DIR);

#get the content of each file to the hash
$out = "$inDirName/appendedClusters.txt";
open(OUT, ">$out") or die "can't open file $out";

foreach my $file (@inFiles){
	if ($file =~ m/^clp/){
		open(IN, "<$inDirName/$file") or die "can't open file $inDirName/$file";
		
		@lines = <IN>;

		while ($lines[0] !~ m/^>/){
			my $a = shift(@lines);
#			print "$a\n";			
		}
		$filesHash{$file} = [@lines];
		
		close(IN);
	}
}

#sort the keys of the hash (file names) by the first line in their content
my @fileNames = keys(%filesHash);

my @sortedFileNames = sort {countOfLine($b) <=> countOfLine($a)} @fileNames;
@sortedFileNames = sort {firstLine($a) cmp firstLine($b)} @sortedFileNames;



#print to the output file
print OUT "------------------------------------------------\n";
foreach my $file (@sortedFileNames){
	print OUT "\t\t\t\t$file:\n\n";
	my @lines = @{$filesHash{$file}};
	print OUT @lines;
	print OUT "\n------------------------------------------------\n";
	
#	print "$file\t$filesHash{$file}->[0]\n";
}

close(OUT);

sub countOfLine{
	my $file = $_[0];
	my $fullLine = uc($filesHash{$file}->[0]);
	if ($fullLine =~ m/>\s*(\d+)\s*X\s*.+/){
		return $1;
	}
}

sub firstLine{
	my $file = $_[0];
	my $fullLine = uc($filesHash{$file}->[0]);
	if ($fullLine =~ m/>\s*\d+\s*X\s*(.+)/){
		return $1;
	}
}