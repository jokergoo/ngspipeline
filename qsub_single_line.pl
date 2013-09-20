#!/bin/perl/bin/

use strict;
use File::Temp qw/tempfile tempdir/;
use File::Basename;

if(scalar(@ARGV) != 2) {
	print <<USAGE;

Send single line command to cluster. The script will first write the
command into a shell script in the temporary directory and then send
it by `qsub`.

Usage:

  perl qsub_single_line.pl "-l walltime=10:00:00,mem=10G" "samtools view ..."

The script expects two arguments, one is setting for qsub and the other is
the single line command. For the settings for qsub, you only need to set 
"-l" and "-N" options. Default value for "-N" is a random string.

Before you run the script, you need to modify the email address and the path of
your temporary directory. Just look into the source code, it is simple.

USAGE
	exit 0;
}

my $qsub_settings = shift(@ARGV);
my $command = join " ", @ARGV;

my $tmp_dir = "/ibios/temp2/guz/clustereo/";
my $email = "z.gu\@dkfz.de";

my ($fh, $filename) = tempfile(DIR => $tmp_dir, SUFFIX => ".sh");

print $fh "#!/bin/sh\n";
my $r = basename($filename);
$r =~s/\.sh$//i;
print $fh "#PBS -N inline_job_$r\n" if($qsub_settings !~/-N\s/);
print $fh "#PBS -j oe\n";
print $fh "#PBS -o $tmp_dir\n";
print $fh "#PBS -M $email\n";
#print $fh "#PBS $qsub_settings\n";

print $fh "\n";
print $fh "$command\n";

close $fh;

print "qsub $qsub_settings $filename\n";
`qsub $qsub_settings $filename`;
