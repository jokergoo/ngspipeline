package CO::NGSPipeline::Pipeline::TopHat;

use strict;
use base qw/CO::NGSPipeline/;

sub new {
	my $class = shift;
	$class = ref($class) ? ref($class) : $class;
	
	my $self = {};
	
	return bless $self, $class;
}

sub run {
	
	my $self = shift;
	
	my $pm = $self->get_pipeline_maker;
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

	###################################################################
	# fastqc
	###################################################################
	$pm->set_job_name("$sample_id"."_tophat_fastqc_r1_$i");
	$qid->{fastqc_r1} = $pipeline->tophat->fastqc(
			fastq      => $r1_fastq,
			output_dir => "$pm->{dir}/fastqc_r1_$i"
	);

	$pm->set_job_name("$sample_id"."_tophat_fastqc_r2_$i");
	$qid->{fastqc_r2} = $pipeline->tophat->fastqc(
			fastq      => $r2_fastq,
			output_dir => "$pm->{dir}/fastqc_r2_$i"
	);

	####################################################################
	# trim
	####################################################################
	$pm->set_job_name("$sample_id"."_tophat_trimmed_$i");
	$qid->{trim} = $pipeline->tophat->trim(
		fastq1  => $r1_fastq,
		fastq2  => $r2_fastq,
		output1 => "$prefix1.trimmed.$i.fastq.gz",
		output2 => "$prefix2.trimmed.$i.fastq.gz",
		polya   => 1,
	);	
				
	###################################################################
	# fastqc after trimming
	###################################################################
	$pm->set_job_name("$sample_id"."_tophat_fastqc_r1_trimmed_$i");
	$pm->set_job_dependency($qid->{trim});
	$qid->{fastqc_r1_trimmed} = $pipeline->tophat->fastqc(
		fastq        => "$prefix1.trimmed.$i.fastq.gz",
		output_dir   => "$pm->{dir}/fastqc_r1_trimmed_$i",
		delete_input => 0
	);

	$pm->set_job_name("$sample_id"."_tophat_fastqc_r2_trimmed_$i");
	$pm->set_job_dependency($qid->{trim});
	$qid->{fastqc_r2_trimmed} = $pipeline->tophat->fastqc(
		fastq        => "$prefix2.trimmed.$i.fastq.gz",
		output_dir   => "$pm->{dir}/fastqc_r2_trimmed_$i",
		delete_input => 0
	);

	##################################################################
	# alignment
	##################################################################
	$pm->set_job_name("$sample_id"."_tophat_align");
	$pm->set_job_dependency($qid->{trim});
	$qid->{align} = $pipeline->tophat->align(
		fastq1 => "$prefix1.trimmed.$i.fastq.gz",
		fastq2 => "$prefix2.trimmed.$i.fastq.gz",
		output => "$pm->{dir}/$sample_id.bam",
		sample_id => $sample_id,
		strand => $is_strand_specific,
		delete_input => 1
	);

	####################################################################
	# more detailed QC
	####################################################################
	$pm->set_job_name("$sample_id"."_rnaseq_qc");
	$pm->set_job_dependency($qid->{align});
	$qid->{qc} = $pipeline->tophat->rnaseqqc(
		bam       => "$pm->{dir}/$sample_id.bam",
		sample_id => $sample_id
	);

	#####################################################################
	# counting and RPKM
	####################################################################
	$pm->set_job_name("$sample_id"."_counting");
	$pm->set_job_dependency($qid->{align});
	$qid->{counting} = $pipeline->tophat->counting(
		bam => "$pm->{dir}/$sample_id.bam",
		strand => $is_strand_specific
	);
}

1;
