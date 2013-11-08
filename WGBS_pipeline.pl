#!/bin/perl/bin/

BEGIN {
	use File::Basename;
	unshift(@INC, dirname($0)."/lib");
}

my $SCRIPT_DIR = dirname($0);

use strict;

use CO::PipelineMaker;
use File::Basename;
use CO::NGSPipeline::Getopt;

use CO::NGSPipeline::Pipeline::Bismark;
use CO::NGSPipeline::Pipeline::BSMAP;
use CO::NGSPipeline::Pipeline::methylCtools;



my $opt = CO::NGSPipeline::Getopt->new;

$opt->before("
WGBS pipeline.

USAGE:

  perl $0 --list file --dir dir --tool tool
  perl $0 --list file --dir dir --tool tool --no-bissnp
  perl $0 --list file --dir dir --tool tool --enforce
  perl $0 --list file --dir dir --tool tool --sample s1,s2
  
");

$opt->after("
NOTE:
  If your fastq files are stored in the standard directory structure which
are generated by data management group, use get_sample_list_from_std_dir.pl
first to generate sample list file.

  Since methylCtools uses BWA to do alignment, about 20% alignment jobs will be
sent to Convey.

FEATURES:
  - record running time for every command
  - catch errors both from exit code and output file size
  - re-run pipeline while skip the upsteam successful jobs
  - generate a detailed QC report

");

my $wd = "analysis";
my $tool = "bsmap";
my $list;
my $std_dir;
my $enforce = 0;
my $request_sampleid;
my $no_bissnp = 0;
my $do_test = 0;
my $filesize = 1024*1024;
my $prefix = "";

$opt->add(\$list, "list=s");
$opt->add(\$wd, "dir=s");
$opt->add(\$tool, "tool=s", "available tools: bismark, bsmap, methyctools");
$opt->add(\$enforce, "enforce");
$opt->add(\$request_sampleid, "sample=s");
$opt->add(\$do_test, "test");
$opt->add(\$filesize, "filesize=i");
$opt->add(\$prefix, "prefix=s");
$opt->add(\$no_bissnp, "no-bissnp", "whether use BisSNP or methylation calling script of each tool to do methylation calling. By default, the three pipelines use BisSNP.");

$opt->getopt;

my $sample = $list;

foreach my $sample_id (sort keys %$sample) {
	
	print "=============================================\n";
	print "submit pipeline for $sample_id\n";
	
	my $r1 = $sample->{$sample_id}->{r1};
	my $r2 = $sample->{$sample_id}->{r2};
	my $library = $sample->{$sample_id}->{library};

	my $pm = CO::PipelineMaker->new(dir => "$wd/$sample_id",
	                                enforce => $enforce,
									do_test => $do_test,
									filesize => $filesize,
									prefix => $prefix,);
	
	my $pipeline;
	if($tool eq "bismark") {
	
		$pipeline = CO::NGSPipeline::Pipeline::Bismark->new();
		
	} elsif($tool eq "bsmap") {
	
		$pipeline = CO::NGSPipeline::Pipeline::BSMAP->new();
		
	} elsif($tool eq "methylctools") {
	
		$pipeline = CO::NGSPipeline::Pipeline::methylCtools->new();
		
	} elsif($tool eq "bsseq") {
	
		$pipeline = CO::NGSPipeline::Pipeline::BSSEQ->new();
		
	} else {
		die "--tool can only be set to one of 'bismark', 'bsmap' and 'methylctools'.\n";
	}
	
	$pipeline->set_pipeline_maker($pm);
	# passing parameters for the pipeline
	$pipeline->run(sample_id => $sample_id,
		               r1 => $r1,
					   r2 => $r2,
					   library => $library,
					   no_bissnp =>$no_bissnp,
					   );
		
}

