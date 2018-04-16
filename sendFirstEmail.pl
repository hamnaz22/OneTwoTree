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

my $from = 'evolseq@post.tau.ac.il';

my $subject = "Your OneTwoTree job $id $jobTitle is being processed.";

my $resultsLink = "http://onetwotree.tau.ac.il/results.html?jobId=$id";
if($jobTitle)
{
	# for the url, replace spaces with %20 in $jobTitle 
	$jobTitle =~ s/ /%20/g;
	$resultsLink .= "&jobTitle=$jobTitle";
}


my $message .= "Hello,\n\n".$subject. "\nA link to your job is available at:\n";
$message .= "$resultsLink\n\n";

$message .= "We will email you again when processing is finished.\n\n";
$message .= "Thank you for using OneTwoTree\n";
#$message .= "OneTwoTree Team\n";

 
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


