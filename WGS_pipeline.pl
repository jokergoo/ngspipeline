#!/bin/perl/bin/

BEGIN {
	use File::Basename;
	unshift(@INC, dirname($0)."/lib");
}

use strict;
use CO::PipelineMaker;
use File::Basename;
use CO::NGSPipeline::Getopt;

# pieline this script supports
use CO::NGSPipeline::Pipeline::WGS;


my $opt = CO::NGSPipeline::Getopt->new;

$opt->before("
ChIPseq pipeline.

USAGE:

  perl $0 --list file --dir dir --tool tool
  perl $0 --list file --dir dir --tool tool --enforce
  perl $0 --list file --dir dir --tool tool --sample s1,s2
  
");

$opt->after("
NOTE:
  If your fastq files are stored in the standard directory structure which
are generated by data management group, use get_sample_list_from_std_dir.pl
first to generate sample list file.

");

# default values
my $wd                 = "analysis";
my $tool               = "wgs";
my $list;
my $enforce            = 0;
my $request_sampleid;
my $do_test            = 0;
my $filesize           = 1024*1024;
my $prefix             = "";
my $email              = 'z.gu@dkfz.de';


$opt->add(\$list,               "list=s");
$opt->add(\$wd,                 "dir=s");
$opt->add(\$tool,               "tool=s", 
                                "available tools: wgs");
$opt->add(\$enforce,            "enforce");
$opt->add(\$request_sampleid,   "sample=s");
$opt->add(\$do_test,            "test");
$opt->add(\$filesize,           "filesize=i");
$opt->add(\$prefix,             "prefix=s");
$opt->add(\$email,              "email=s");

$opt->getopt;



foreach my $sample_id (sort keys %$list) {
	
	my $r1 = $list->{$sample_id}->{r1};
	my $r2 = $list->{$sample_id}->{r2};
	
	#if(($tool eq "defuse" or $tool eq "fusionmap" or $tool eq "fusionhunter") and scalar(@$r1) > 1) {
	#	die "Currently only support one lane per sample for gene fusion pipeline. $sample_id has multiple lanes.\n";
	#}
}

foreach my $sample_id (sort keys %$list) {
	
	print "=============================================\n";
	print "submit pipeline for $sample_id\n";
	
	my $pipeline;
	if($tool eq "wgs") {
	
		$pipeline = CO::NGSPipeline::Pipeline::WGS->new();
		
	} else {
		die "--tool can only be set to one of 'wgs'.\n";
	}
	
	my $r1 = $list->{$sample_id}->{r1};
	my $r2 = $list->{$sample_id}->{r2};
	my $library = $list->{$sample_id}->{library};

	my $pm = CO::PipelineMaker->new(dir      => "$wd/$sample_id",
	                                enforce  => $enforce,
									do_test  => $do_test,
									filesize => $filesize,
									prefix   => $prefix,
									email    => $email,);

	$pipeline->set_pipeline_maker($pm);
	$pipeline->run(sample_id              => $sample_id,
		               r1                 => $r1,
					   r2                 => $r2,
					   library            => $library,
					   );
					   
}

