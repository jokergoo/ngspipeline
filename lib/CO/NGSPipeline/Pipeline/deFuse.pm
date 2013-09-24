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
	
	# prefix means absolute path without fast/fq or fast.gz/fq.gz
	my $prefix1 = basename($r1->[0]);
	$prefix1 =~s/\.(fq|fastq)(\.gz)?$//;
	$prefix1 = "$pm->{dir}/$prefix1";
	my $prefix2 = basename($r2->[0]);
	$prefix2 =~s/\.(fq|fastq)(\.gz)?$//;
	$prefix2 = "$pm->{dir}/$prefix2";
	
	
	
	my $qid;

	my $r1_fastq = $r1->[0];
	my $r2_fastq = $r2->[0];
	
	
	####################################################################
	# trim
	####################################################################
	$pm->set_job_name("$sample_id"."_defuse_trimmed");
	$qid->{trim} = $pipeline->defuse->trim(
		fastq1  => $r1_fastq,
		fastq2  => $r2_fastq,
		output1 => "$prefix1.trimmed.fastq.gz",
		output2 => "$prefix2.trimmed.fastq.gz",
		polya   => 1,
	);	
				
	$pm->set_job_name("$sample_id"."_defuse");
	$qid = $pipeline->genefusion->defuse(
		fastq1 => "$prefix1.trimmed.fastq.gz",
		fastq2 => "$prefix2.trimmed.fastq.gz",
		delete_input => 1,
	);
}

1;
