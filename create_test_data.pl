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
Extract n sequences from the FASTQ file.

Usage:

  perl create_testdata.pl --input fastq --output output --start start --records n
  
  --input   FastQ file
  --output  test FastQ file
  --start   Which read to start
  --records How many reads do you want
  
USAGE
	
	exit 0;
}

use CO::FastQ;

my $fastq;
my $output;
my $start = 1;
my $n_records = 10000;

GetOptions( "input=s" => \$fastq,
            "output=s" => \$output;
			"start=i" => \$start,
			"records=i" => \$n_records) or die;
			
my $fastq = CO::FastQ->new(file => $fastq);

my $fh;
if(-p $output) {
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

