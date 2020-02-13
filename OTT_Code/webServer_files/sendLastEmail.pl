use strict;
use warnings;

use Getopt::Long;

my $jobTitle;
my $toEmail;
my $id;

GetOptions(
	"jobTitle=s"        => \$jobTitle,
	"toEmail=s"        => \$toEmail,
	"id=s"     => \$id
);

my $from = 'evolseq@tauex.tau.ac.il';

my $subject;
my $subject_fail = "Your OneTwoTree job $id $jobTitle has failed";
my $final_status_file = "/bioseq/data/results/oneTwoTree/$id/SummaryDir/FinalStatus.txt";
my $final_pass_file = "/bioseq/data/results/oneTwoTree/$id/JOB_PASS.txt";

my $resultsLink = "http://onetwotree.tau.ac.il/results.html?jobId=$id";
my $message .= "Hello,\n\n";
if (-e $final_pass_file)
{
	$subject = "Your OneTwoTree results of run $id $jobTitle are ready";
	$message .= "$subject at:\n";
	$message .= "Your job was completed successfully\n";
} else {
	$subject = "Your OneTwoTree run $id $jobTitle has failed";
	$message .= "$subject_fail at:\n";
	$message .= "Please contact us for more info\n";
}


if($jobTitle)
{
	# for the url, replace spaces with %20 in $jobTitle 
	$jobTitle =~ s/ /%20/g;
	$resultsLink .= "&jobTitle=$jobTitle";
}



#my $message .= "Hello,\n\n".$subject. " at:\n";
$message .= "$resultsLink\n\n";


$message .= "Please note: the results will be kept on the server for one month.\n\n";
$message .= "Thank you for using OneTwoTree,\n";
$message .= "OneTwoTree Team\n";

 
#open(MAIL, "|/usr/sbin/sendmail -t");
open(MAIL, "|/usr/sbin/sendmail -t") or die "can't open mail";
 
# Email Header
print MAIL "To: $toEmail\n";
print MAIL "From: $from\n";
print MAIL "Subject: $subject\n\n";

# Email Body
print MAIL $message;
my $return_code = $? >> 8;
print "print MAIL returned: $return_code\n";

close(MAIL);
#print "Email Sent Successfully to: $toEmail\n";

#use Sys::Hostname;
#$host = hostname;
my $host = `hostname`;

print "hostname: $host\n";


