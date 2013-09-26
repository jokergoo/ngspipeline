#!/bin/perl/bin
use strict;

BEGIN {
	use File::Basename;
	unshift(@INC, dirname($0)."/lib");
}

###########################
# usage:
#   perl create_testdata.pl --input fastq --output output --start start --records n

if(! scalar(@ARGV)) {
	print <<USAGE;
Extract n sequences from FastQ file.

Usage:

  perl create_testdata.pl -i fastq -o output -s start -k n
  
  --input, -i   FastQ file, gzipped or not.
  --output, -o  output FastQ file. gzipped or not (default: uncompressed FastQ
                to STDOUT)
  --start, -s   Which read to start (default: 1)
  --records, -k How many reads do you want (default: 10000)
  
USAGE
	
	exit 0;
}

use CO::FastQ;
use Getopt::Long;
use Fcntl;

my $fastq;
my $output;
my $start = 1;
my $n_records = 10000;

GetOptions( "input|i=s" => \$fastq,
            "output|o=s" => \$output,
			"start|s=i" => \$start,
			"records|k=i" => \$n_records) or die;
			
my $fastq = CO::FastQ->new(file => $fastq);

my $fh;
if((!defined($output)) or $output eq "-") {
	$fh = \*STDOUT;
} elsif(-p $output) {
	sysopen($fh, "$output", O_WRONLY);
} elsif($output =~/\.gz$/) {
	open $fh, " | gzip -c > $output" or die "cannot create pipe to gzip\n";
} else {
	open $fh, ">$output" or die "cannot create $output\n";
}
	
my $i_processed = 0;
while(my $read = $fastq->next) {
	
	if($fastq->i < $start) {
		next;
	}
	
	print $fh $read->record;
	
	if($i_processed >= $n_records) {
		last;
	}
}
close $fh;

