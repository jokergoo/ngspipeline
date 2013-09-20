#!/bin/perl/bin
use strict;

###########################
# usage:
#   perl createtestdata.pl -i fastq -s start -k n_records -o output

if(! scalar(@ARGV)) {
	print <<USAGE;
Extract the first n sequences from the FASTQ file.

Usage:

  zcat xx.fq.gz | perl create_test_data.pl -n 100000 | gzip > test.fq.gz
  
  -n Number of sequences to extract, default is 20
  -f FASTQ file. If not -f option is set, it will read from STDIN.

The script expects uncompressed fastq file and write fastq file to STDOUT, 
so you should put the command in pipes.
  
USAGE
	
	exit 0;
}


