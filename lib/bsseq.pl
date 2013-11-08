

use strict;

use CO::PipelineMaker;
use File::Basename;

use CO::NGSPipeline::Pipeline::BSSEQ;


open F, "list_bedgraph";

my $sample;
while(my $line = <F>) {
	chomp $line;
	my @tmp = split "\t", $line;
	$sample->{$tmp[0]} = $tmp[1];
}


foreach my $sample_id (sort keys %$sample) {
	
	print "=============================================\n";
	print "submit pipeline for $sample_id\n";
	
	mkdir "/icgc/lsdf/mb/analysis/B060_016/WGBS_analysis/analysis_bsmap_combine_dmr/$sample_id" unless(-e "/icgc/lsdf/mb/analysis/B060_016/WGBS_analysis/analysis_bsmap_combine_dmr/$sample_id");

	my $pm = CO::PipelineMaker->new(dir => "/icgc/lsdf/mb/analysis/B060_016/WGBS_analysis/analysis_bsmap_combine_dmr/$sample_id",
	                                enforce => 1,
									do_test => 0,
									filesize => 1024*0124,
									prefix => "",);
	
	my $pipeline = CO::NGSPipeline::Pipeline::BSSEQ->new();
		
	$pipeline->set_pipeline_maker($pm);
	# passing parameters for the pipeline
	$pipeline->run(sample_id => $sample_id,
		            bedgraph => $sample->{$sample_id});
		
}

