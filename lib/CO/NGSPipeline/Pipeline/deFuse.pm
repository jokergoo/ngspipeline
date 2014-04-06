package CO::NGSPipeline::Pipeline::deFuse;

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
	             r1 => undef,
				 r2 => undef,
				 library => undef,
				 dir => undef,
				 is_strand_specific => 0,
				 @_);
				 
	my $sample_id = $param{sample_id};
	my $r1 = $param{r1};
	my $r2 = $param{r2};
	my $library = $param{library};
	my $dir = $param{dir};
	my $is_strand_specific = $param{is_strand_specific};
	
		
	my $qid;
	my $trimmed_fastq1_arrayref = [];
	my $trimmed_fastq2_arrayref = [];
	
	for(my $i = 0; $i < scalar(@$r1); $i ++) {
		
		my $r1_fastq = $r1->[$i];
		my $r2_fastq = $r2->[$i];
		
		
		####################################################################
		# trim
		####################################################################
		$pm->set_job_name("$sample_id"."_defuse_trimmed");
		$qid->{trim}->[$i] = $pipeline->genefusion->trim(
			fastq1  => $r1_fastq,
			fastq2  => $r2_fastq,
			output1 => "$pm->{dir}/$sample_id.$i.trimmed.r1.fastq.gz",
			output2 => "$pm->{dir}/$sample_id.$i.trimmed.r2.fastq.gz",
			polya   => 1,
		);
		
		push(@$trimmed_fastq1_arrayref, "$pm->{dir}/$sample_id.$i.trimmed.r1.fastq.gz");
		push(@$trimmed_fastq2_arrayref, "$pm->{dir}/$sample_id.$i.trimmed.r2.fastq.gz");
		
	}			
		
	$pm->set_job_name("$sample_id"."_defuse");
	$pm->set_job_dependency(@{$qid->{trim}});
	$qid = $pipeline->genefusion->defuse(
		fastq1_arrayref => $trimmed_fastq1_arrayref,
		fastq2_arrayref => $trimmed_fastq2_arrayref,
		delete_input => 1,
	);

}

1;
