package CO::NGSPipeline::Pipeline::STAR;

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

	print $sample_id;

	###################################################################
	# fastqc
	###################################################################
	$pm->set_job_name("$sample_id"."_star_fastqc_r1");
	$qid->{fastqc_r1} = $pipeline->star->fastqc(
			fastq      => $r1_fastq,
			output_dir => "$pm->{dir}/fastqc_r1"
	);

	$pm->set_job_name("$sample_id"."_star_fastqc_r2");
	$qid->{fastqc_r2} = $pipeline->star->fastqc(
			fastq      => $r2_fastq,
			output_dir => "$pm->{dir}/fastqc_r2"
	);

	####################################################################
	# trim
	####################################################################
	$pm->set_job_name("$sample_id"."_star_trimmed");
	$qid->{trim} = $pipeline->star->trim(
		fastq1  => $r1_fastq,
		fastq2  => $r2_fastq,
		output1 => "$prefix1.trimmed.fastq.gz",
		output2 => "$prefix2.trimmed.fastq.gz",
		polya   => 1,
	);	
				
	###################################################################
	# fastqc after trimming
	###################################################################
	$pm->set_job_name("$sample_id"."_star_fastqc_r1_trimmed");
	$pm->set_job_dependency($qid->{trim});
	$qid->{fastqc_r1_trimmed} = $pipeline->star->fastqc(
		fastq        => "$prefix1.trimmed.fastq.gz",
		output_dir   => "$pm->{dir}/fastqc_r1_trimmed",
		delete_input => 0
	);

	$pm->set_job_name("$sample_id"."_star_fastqc_r2_trimmed");
	$pm->set_job_dependency($qid->{trim});
	$qid->{fastqc_r2_trimmed} = $pipeline->star->fastqc(
		fastq        => "$prefix2.trimmed.fastq.gz",
		output_dir   => "$pm->{dir}/fastqc_r2_trimmed",
		delete_input => 0
	);

	##################################################################
	# alignment
	##################################################################
	$pm->set_job_name("$sample_id"."_star_align");
	$pm->set_job_dependency($qid->{trim});
	$qid->{align} = $pipeline->star->align(
		fastq1       => "$prefix1.trimmed.fastq.gz",
		fastq2       => "$prefix2.trimmed.fastq.gz",
		output       => "$pm->{dir}/$sample_id.bam",
		sample_id    => $sample_id,
		strand       => $is_strand_specific,
		delete_input => 1,
	);
	
	####################################################################
	# flagstat
	####################################################################
	$pm->set_job_name("$sample_id"."_star_flagstat");
	$pm->set_job_dependency($qid->{align});
	$qid->{flagstat} = $pipeline->star->samtools_flagstat(
		sam          => "$pm->{dir}/$sample_id.bam",
		output       => "$pm->{dir}/sample_id.flagstat",
		delete_input => 0
	);

	####################################################################
	# more detailed QC
	####################################################################
	$pm->set_job_name("$sample_id"."_rnaseq_qc");
	$pm->set_job_dependency($qid->{align});
	$qid->{qc} = $pipeline->star->rnaseqqc(
		bam       => "$pm->{dir}/$sample_id.bam",
		sample_id => $sample_id
	);

	#####################################################################
	# counting and RPKM
	####################################################################
	$pm->set_job_name("$sample_id"."_counting");
	$pm->set_job_dependency($qid->{align});
	$qid->{counting} = $pipeline->star->counting(
		bam => "$pm->{dir}/$sample_id.bam",
		strand => $is_strand_specific
	);
}

1;
