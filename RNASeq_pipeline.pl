#!/bin/perl/bin/

BEGIN {
	use File::Basename;
	unshift(@INC, dirname($0)."/lib");
}

my $SCRIPT_DIR = dirname($0);

use strict;
use Getopt::Long;
use CO::NGSPipeline::Utils;
use CO::NGSPipeline;
use File::Basename;
use List::Vectorize;

if(! scalar(@ARGV)) {
	print <<USAGE;
########################################################################
#  Please confirm you have changed email in lib/NGSPipeline/Config.pm  #
########################################################################

USAGE:
  
  perl $0 --list file --dir dir --tool tool
  perl $0 --list file --dir dir --tool tool --enforce
  perl $0 --list file --dir dir --tool tool --sample s1,s2
  
  --list      sample list, containing columns which are:
              1. fastq file for paired end 1, should be gzipped
              2. fastq file for paired end 2, should be gzipped
              3. sample name. unique

  --dir       working dir, default is `analysis`. Under the working dir, there
              are list of directories named with sample names which are called
              job directory for each sample.
  
  --tool      in "gsnap", "tophat", "star"
              "defuse", "tophatfusion", "fusionmap", "fusionhunter"
      
  --enforce   enforce to re-run pipeline from the beginning no matter they were
              successfully finished or not.
  
  --sample    subset of sample ids, should seperated by "," (no blank)
		
  --filesize  If size of some output files (e.g. bam files, methylation calling 
              files) are smaller than this value, then step is terminated. 
              Default is 1M (1024*1024). Set it to 0 or non-number
              strings to shut down file size checking.
  --strand    strand specific
  	
	exit 0;
}

my $wd = "analysis";
my $tool = "";
my $list;
my $std_dir;
my $enforce = 0;
my $request_sampleid;
my $do_test = 0;
my $filesize = 1024*1024;
my $is_strand_specific = 0;

GetOptions("list=s"    => \$list,
           "dir=s"     => \$wd,
           "tool=s"    => \$tool,
		   "enforce"   => \$enforce,
		   "sample=s"  => \$request_sampleid,
		   "filesize"  => \$filesize,
		   "test"      => \$do_test,
		   "strand"    => \$is_strand_specific);

		   
my %subset_samples = map { $_ => 1} split ",", $request_sampleid;
$filesize += 0;

open F, $list or die "Cannot open $list\n";
my $r1;
my $r2;
my $sample;
my $n_sample = 0;
while(my $line = <F>) {
	chomp $line;
	next if($line =~/^\s*$/);
	next if($line =~/^#/);
	
	my @tmp = split "\t", $line;
	$tmp[0] = to_abs_path($tmp[0]);
	$tmp[1] = to_abs_path($tmp[1]);
	
	if(basename($tmp[0]) eq basename($tmp[1])) {
		die "two fastq files have same names! check your file\n";
	}
	
	if(scalar(%subset_samples) and !$subset_samples{$tmp[2]}) {
		print "$tmp[2] is not in --sample, skip this sample.\n";
		next;
	}
	
	# if no record for this sample, initialize the array reference
	if(! defined($sample->{$tmp[2]})) {
		$sample->{$tmp[2]} = {};
		$sample->{$tmp[2]}->{r1} = [];
		$sample->{$tmp[2]}->{r2} = [];
		$sample->{$tmp[2]}->{library} = [];
	}
	
	push(@{$sample->{$tmp[2]}->{r1}}, $tmp[0]);
	push(@{$sample->{$tmp[2]}->{r2}}, $tmp[1]);
	
	# currently do not support multiple libraries for a same sample
	push(@{$sample->{$tmp[2]}->{library}}, defined($tmp[3]) ? $tmp[3] : undef);
	
	$n_sample ++;
}

$wd = to_abs_path($wd);
# seems set mode of the dir to 0755 not always successful
-e $wd ? 1: mkdir $wd, 0775 || die "cannto create dir: $wd with mode 0775\n";

$tool = lc($tool);

print "Working directory is $wd.\n";
print "Totally $n_sample samples with ". scalar(keys %$sample)." unique sample ids \n";
print "Using $tool.\n\n";

my $command = [];
foreach my $sample_id (sort keys %$sample) {
	
	print "=============================================\n";
	print "submit pipeline for $sample_id\n";
	
	my $r1 = $sample->{$sample_id}->{r1};
	my $r2 = $sample->{$sample_id}->{r2};
	my $library = $sample->{$sample_id}->{library};

	my $pipeline = NGSPipeline->new(dir => "$wd/$sample_id",
	                                enforce => $enforce,
									do_test => $do_test,
									filesize => $filesize);

	# prefix means absolute path without fast/fq or fast.gz/fq.gz
	my $prefix1 = basename($r1->[0]);
	$prefix1 =~s/\.(fq|fastq)(\.gz)?$//;
	$prefix1 = "$pipeline->{dir}/$prefix1";
	my $prefix2 = basename($r2->[0]);
	$prefix2 =~s/\.(fq|fastq)(\.gz)?$//;
	$prefix2 = "$pipeline->{dir}/$prefix2";
	
	if($tool eq "gsnap") {
		
		require("$SCRIPT_DIR/pipeline/gsnap.pl");
		
	} elsif($tool eq "tophat") {
		
		require("$SCRIPT_DIR/pipeline/tophat.pl");
		
	} elsif($tool eq "star") {
		
		require("$SCRIPT_DIR/pipeline/star.pl");
										 
	} elsif($tool eq "fusionmap") {
	
		require("$SCRIPT_DIR/pipeline/fusionmap.pl");
		
	} elsif($tool eq "defuse") {
		
		require("$SCRIPT_DIR/pipeline/defuse.pl");
		
	} elsif($tool eq "tophatfusion" or $tool eq "tophat-fusion") {
		
		require("$SCRIPT_DIR/pipeline/tophatfusion.pl");
		
	} elsif($tool eq "fusionhunter") {
		
		require("$SCRIPT_DIR/pipeline/fusionhunter.pl");
		
	} else {
		die "--tool can only be set to one of 'gsnap', 'tophat', 'star', 'fusionmap', 'defuse', 'fusionhunter'.\n";
	}
}

