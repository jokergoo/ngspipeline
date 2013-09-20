#!/bin/perl/bin/

BEGIN {
	use File::Basename;
	unshift(@INC, dirname($0)."/../lib");
}

use Template;
use File::Temp qw/tempfile tempdir/;
use File::Find;
use strict;
use List::Vectorize;
use File::Basename;
use Getopt::Long;
use Data::Dumper;


if(! scalar(@ARGV)) {
	print <<USAGE;
Usage:
	perl bsmap_report.pl --tool tool --dir /dir/of/the/sample --sample sample_name
	
    If --sample is not specified, it is basename of the value of --dir

USAGE
exit;
}

my $sample_dir;
my $sample_id;
my $tool;

GetOptions("tool=s" => \$tool,
           "dir=s" => \$sample_dir,
           "sample=s" => \$sample_id);

unless($tool eq "bismark" or $tool eq "bsmap" or $tool eq "methylctools") {
	die "--tool should only be in 'bismark', 'bsmap' and 'methylctools'\n";
}
      
if(! $sample_id) {
	$sample_id = basename($sample_dir);
}
		   
my $script_dir = dirname($0);

my $tt = Template->new({
	INCLUDE_PATH => ["$script_dir/tt/src",
	                 "$script_dir/tt/lib"]
	});

my $tmp_dir = "$sample_dir/tmp";

# find how many fastqc files
my @fastqc_dir = grep { /\/fastqc_r\d+(_trimmed)?_\d+$/} glob("$sample_dir/*");
our $fastqc_data = [];
find(sub {/fastqc_data\.txt$/ and push(@$fastqc_data, $File::Find::name)}, @fastqc_dir);

# insert size files
my ($insertsize_file) = glob("$sample_dir/*.insertsiz.metrics");

# flagstat files
my $flagstat_file = [glob("$sample_dir/*.flagstat")];  # maybe more than one flagstat file

# mark duplicate files
my $duplicate_file = [glob("$sample_dir/*.mkdup.metrics")];

########################################################
# read cpg.bedgraph from bissnp
my ($bissnp_cpg_bedgraph_file) = glob("$sample_dir/*.CG.bedgraph");

my ($lambda_conversion_rate_file) = glob("$sample_dir/*.lambda.conversion.txt");

my $input = 'report_bs.tt';
my $vars  = {
    sample_id => $sample_id,
    fastqc_file => $fastqc_data,
    insertsize_file => $insertsize_file,
    flagstat_file => $flagstat_file,
	duplicate_file => $duplicate_file,
	bissnp_cpg_bedgraph_file => $bissnp_cpg_bedgraph_file,
	lambda_conversion_rate_file => $lambda_conversion_rate_file,
	tool => $tool,
};

#print Dumper $vars;

use Template::Stash;

# define list method to return new list of odd numbers only
$Template::Stash::SCALAR_OPS->{ "rand" } = sub {
    return int(rand(1000000));
};
$Template::Stash::SCALAR_OPS->{ "basename" } = sub {
    return basename($_[0]);
};


my ($fh, $filename) = tempfile(DIR => $tmp_dir, SUFFIX => ".Rnw");
close $fh;
print "tmp Rsweave file: $filename\n";
$tt->process($input, $vars, $filename) || die $tt->error( );

system("cd $tmp_dir;R CMD Sweave --pdf $filename");
$filename =~s/\.Rnw/.pdf/;
system("cp $filename $sample_dir/$sample_id.bsmap.report.pdf");

print "report: $sample_dir/$sample_id.bsmap.report.pdf\n";


