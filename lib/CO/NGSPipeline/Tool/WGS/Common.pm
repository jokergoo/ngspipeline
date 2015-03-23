package CO::NGSPipeline::Tool::WGS::Common;

##############################################################################
# provide command method for each pipeline. It is a base module for all specific
# pipelines.

use strict;
use CO::NGSPipeline::Tool::WGS::Config;
use CO::NGSPipeline::Tool::Config;
use CO::Utils;
use File::Basename;

use base qw(CO::NGSPipeline::Tool::Common
			CO::NGSPipeline);

sub new {
	my $class = shift;
	$class = ref($class) ? ref($class) : $class;

	my $self = {};
	
	return bless $self, $class;
}

sub coverage {
	my $self = shift;
	
	my %param = ( "sam" => undef,
				"ungap" => undef,
				  "output" => undef,
				  @_);
	
	my $sam    = to_abs_path( $param{sam} );
	my $sample_id = $param{sample_id};
	my $ungap = to_abs_path( $param{ungap} );
	my $output = to_abs_path( $param{output} );
	my $delete_input = $param{delete_input};
	
	my $dir = dirname($output);

	my $pm = $self->get_pipeline_maker;
	
	$pm->add_command("/home/guz/co_pipeline/Roddy/Plugins/COWorkflows/resources/analysisTools/qcPipelineTools/coverageQcD/coverageQc --alignmentFile=$sam --outputFile=$output --processors=4 --basequalCutoff=25 --ungappedSizes=$ungap");
	$pm->add_command("python /home/guz/co_pipeline/Roddy/Plugins/COWorkflows/resources/analysisTools/qcPipeline/genomeCoverage.py --alignmentFile=$sam,$dir/$sample_id.readCoverage_1kb_windows.txt --processors=12 --countReads --windowSize=1 --ignore_chrRandom_chrM_hap_chrUn");

	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_common_coverage",
							 "-l" => { nodes => "1:ppn=12:lsdf", 
									    mem => "2GB",
										walltime => "50:00:00"});
	return($qid);
}

1;
