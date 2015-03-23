package CO::NGSPipeline::Tool::BSseq::Common;


use strict;
use CO::NGSPipeline::Tool::Config;
use CO::Utils;
use File::Basename;
use File::Temp qw/tempfile/;

use File::Basename;
use CO::Utils;
use base qw/CO::NGSPipeline/;

sub bsqc {
	my $self = shift;
	
	my %param = ( "dir" => undef,
	              "tool" => undef,
	              "sample" => undef,
				  @_);
	
	my $sample_dir    = to_abs_path( $param{dir} );
	my $tool = $param{tool};
	my $sample = $param{sample};
	my $base_dir = to_abs_path( $param{base_dir});
	
	my $pm = $self->get_pipeline_maker;
	
	my $script_dir = to_abs_path(dirname(__FILE__)."/../../Report");
	$pm->add_command("perl $script_dir/bs_report.pl --tool $tool --dir $sample_dir --sample $sample");

	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_common_bsqc",
							 "-l" => { nodes => "1:ppn=1:lsdf", 
									    mem => "20GB",
										walltime => "20:00:00"});
	return($qid);

}



1;
