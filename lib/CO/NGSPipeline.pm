package CO::NGSPipeline;

use strict;

# This module is a universal base module for all pipeline.
# Its main function is to pass the PipelineMaker object to each pipeline and 
# each step

sub set_pipeline_maker {
	my $pipeline = shift;
	my $pm = shift;
	
	$pipeline->{_pipeline_maker} = $pm;
	
	return $pipeline;
}

sub get_pipeline_maker {
	my $pipeline = shift;
	
	return $pipeline->{_pipeline_maker};
}

1;
