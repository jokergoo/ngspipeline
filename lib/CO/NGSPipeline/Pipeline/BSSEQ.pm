package CO::NGSPipeline::Pipeline::BSSEQ;

use strict;
use File::Basename;
use base qw/CO::NGSPipeline/;

sub new {
	my $class = shift;
	$class = ref($class) ? ref($class) : $class;
	
	my $self = {};
	
	return bless $self, $class;
}

sub run {
	
	my $pipeline = shift;
	
	my $pm = $pipeline->get_pipeline_maker;
	unless(defined($pm)) {
		die "pipeline maker should be defined in a pipeline.\n";
	}
	
	my %param = (sample_id => undef,
	             bedgraph => undef,
	             @_);
				 
	my $sample_id = $param{sample_id};
	my $bedgraph = $param{bedgraph};
	
	my $qid = {};
	for my $chr (map {"chr$_"} (1..22, "X", "Y")) {
		$pm->set_job_name("$sample_id"."_bsseq_RData_$chr");
		$qid->{RData} = $pipeline->bsmooth->RData(
			input   => "$bedgraph",
			sample_id  => $sample_id,
			chr        => $chr,
			output     => "$pm->{dir}/bsseq2_$chr.RData",
		);
	}
}

1;
