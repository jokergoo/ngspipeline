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
				 remove_duplicate => 0,
				 @_);
				 
	my $sample_id = $param{sample_id};
	my $r1 = $param{r1};
	my $r2 = $param{r2};
	my $library = $param{library};
	my $dir = $param{dir};
	my $is_strand_specific = $param{is_strand_specific};
	my $remove_duplicate = $param{remove_duplicate};
	
	# prefix means absolute path without fast/fq or fast.gz/fq.gz
	my $prefix1 = basename($r1->[0]);
	$prefix1 =~s/\.(fq|fastq)(\.gz)?$//;
	$prefix1 = "$pm->{dir}/$prefix1";
	my $prefix2 = basename($r2->[0]);
	$prefix2 =~s/\.(fq|fastq)(\.gz)?$//;
	$prefix2 = "$pm->{dir}/$prefix2";
	
	

	my $qid = {};
	$qid->{alignment} = [];
	my $sam_sort_list = [];
	for(my $i = 0; $i < scalar(@$r1); $i ++) {

		my $r1_fastq = $r1->[$i];
		my $r2_fastq = $r2->[$i];

		###################################################################
		# fastqc
		###################################################################
		$pm->set_job_name("$sample_id"."_star_fastqc_r1_$i");
		$qid->{fastqc_r1} = $pipeline->star->fastqc(
				fastq      => $r1_fastq,
				output_dir => "$pm->{dir}/fastqc_r1_$i"
		);

		$pm->set_job_name("$sample_id"."_star_fastqc_r2_$i");
		$qid->{fastqc_r2} = $pipeline->star->fastqc(
				fastq      => $r2_fastq,
				output_dir => "$pm->{dir}/fastqc_r2_$i"
		);

=pod
		####################################################################
		# trim
		####################################################################
		$pm->set_job_name("$sample_id"."_star_trimmed_$i");
		$qid->{trim} = $pipeline->star->trim(
			fastq1  => $r1_fastq,
			fastq2  => $r2_fastq,
			output1 => "$prefix1.trimmed.$i.fastq.gz",
			output2 => "$prefix2.trimmed.$i.fastq.gz",
		);	
					
		###################################################################
		# fastqc after trimming
		###################################################################
		$pm->set_job_name("$sample_id"."_star_fastqc_r1_trimmed_$i");
		$pm->set_job_dependency($qid->{trim});
		$qid->{fastqc_r1_trimmed} = $pipeline->star->fastqc(
			fastq        => "$prefix1.trimmed.$i.fastq.gz",
			output_dir   => "$pm->{dir}/fastqc_r1_trimmed_$i",
			delete_input => 0
		);

		$pm->set_job_name("$sample_id"."_star_fastqc_r2_trimmed");
		$pm->set_job_dependency($qid->{trim});
		$qid->{fastqc_r2_trimmed} = $pipeline->star->fastqc(
			fastq        => "$prefix2.trimmed.fastq.gz",
			output_dir   => "$pm->{dir}/fastqc_r2_trimmed",
			delete_input => 0
		);
=cut
		##################################################################
		# alignment
		##################################################################
		$pm->set_job_name("$sample_id"."_star_align_$i");
		$pm->set_job_dependency($qid->{trim});
		$qid->{alignment}->[$i] = $pipeline->star->align(
			fastq1       => "$r1_fastq",
			fastq2       => "$r2_fastq",
			output       => "$pm->{dir}/$sample_id.$i.bam",
			sample_id    => $sample_id,
			strand       => $is_strand_specific,
			delete_input => 0,
		);
		
		####################################################################
		# flagstat
		####################################################################
		$pm->set_job_name("$sample_id"."_star_flagstat_$i");
		$pm->set_job_dependency($qid->{alignment}->[$i]);
		$qid->{flagstat} = $pipeline->star->samtools_flagstat(
			sam          => "$pm->{dir}/$sample_id.$i.bam",
			output       => "$pm->{dir}/$sample_id.$i.flagstat",
			delete_input => 0
		);
		
		$sam_sort_list->[$i] = "$pm->{dir}/$sample_id.$i.bam";
	}
	
	if($remove_duplicate) {
		$pm->set_job_name("$sample_id"."_star_merge_and_nodup");
		$pm->set_job_dependency(@{$qid->{alignment}});
		$qid->{remove_duplicate} = $pipeline->star->merge_nodup(
			sam_list     => $sam_sort_list,
			output       => "$pm->{dir}/$sample_id.merged.bam",
			delete_input => 1
		);
	} else {
		$pm->set_job_name("$sample_id"."_star_merge");
		$pm->set_job_dependency(@{$qid->{alignment}});
		$qid->{remove_duplicate} = $pipeline->star->merge_sam(
			sam_list     => $sam_sort_list,
			output       => "$pm->{dir}/$sample_id.merged.bam",
			delete_input => 1
		);
	}
	
	####################################################################
	# more detailed QC
	####################################################################
	$pm->set_job_name("$sample_id"."_rnaseq_qc");
	$pm->set_job_dependency($qid->{remove_duplicate});
	$qid->{qc} = $pipeline->star->rnaseqqc(
		bam       => "$pm->{dir}/$sample_id.merged.bam",
		sample_id => $sample_id
	);

	#####################################################################
	# counting and RPKM
	####################################################################
	$pm->set_job_name("$sample_id"."_counting");
	$pm->set_job_dependency($qid->{remove_duplicate});
	$qid->{counting} = $pipeline->star->counting(
		bam => "$pm->{dir}/$sample_id.merged.bam",
		strand => $is_strand_specific
	);
}

1;
