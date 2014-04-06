package CO::NGSPipeline::Pipeline::TopHat;

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
	
	my %param = (sample_id          => undef,
	             r1                 => undef,
				 r2                 => undef,
				 library            => undef,
				 dir                => undef,
				 is_strand_specific => 0,
				 remove_duplicate   => 0,
				 @_);
				 
	my $sample_id          = $param{sample_id};
	my $r1                 = $param{r1};
	my $r2                 = $param{r2};
	my $library            = $param{library};
	my $dir                = $param{dir};
	my $is_strand_specific = $param{is_strand_specific};
	my $remove_duplicate   = $param{remove_duplicate};


	my $qid = {};
	$qid->{flagstat} = [];
	my $sam_sort_list = [];
	for(my $i = 0; $i < scalar(@$r1); $i ++) {

		my $r1_fastq = $r1->[$i];
		my $r2_fastq = $r2->[$i];

		###################################################################
		# fastqc
		###################################################################
		#$pm->set_job_name("$sample_id"."_tophat_fastqc_r1_$i");
		#$qid->{fastqc_r1} = $pipeline->tophat->fastqc(
	#			fastq      => $r1_fastq,
	#			output_dir => "$pm->{dir}/fastqc_r1_$i"
	#	);

	#	$pm->set_job_name("$sample_id"."_tophat_fastqc_r2_$i");
	#	$qid->{fastqc_r2} = $pipeline->tophat->fastqc(
	#			fastq      => $r2_fastq,
	#			output_dir => "$pm->{dir}/fastqc_r2_$i"
	#	);

		##################################################################
		# alignment
		##################################################################
		$pm->set_job_name("$sample_id"."_tophat_align_$i");
		$pm->set_job_dependency($qid->{trim});
		$qid->{alignment} = $pipeline->tophat->align(
			fastq1       => "$r1_fastq",
			fastq2       => "$r2_fastq",
			output       => "$pm->{dir}/$sample_id.$i.bam",
			sample_id    => $sample_id,
			strand       => $is_strand_specific,
			delete_input => 0,  # don't delete original fastq files
		);
		
		####################################################################
		# flagstat
		####################################################################
		$pm->set_job_name("$sample_id"."_tophat_flagstat_$i");
		$pm->set_job_dependency($qid->{alignment});
		$qid->{flagstat}->[$i] = $pipeline->tophat->samtools_flagstat(
			sam          => "$pm->{dir}/$sample_id.$i.bam",
			output       => "$pm->{dir}/$sample_id.$i.flagstat",
			delete_input => 0,
		);
		
		$sam_sort_list->[$i] = "$pm->{dir}/$sample_id.$i.bam";
	}
	
	$pm->set_job_name("$sample_id"."_tophat_merge_and_mkdup");
	$pm->set_job_dependency(@{$qid->{flagstat}});
	$qid->{remove_duplicate} = $pipeline->tophat->merge_nodup(
			sam_list     => $sam_sort_list,
			output       => "$pm->{dir}/$sample_id.mkdup.bam",
			delete_input => 1,
			REMOVE_DUPLICATES => 'FALSE',
	);

	$pm->set_job_name("$sample_id"."_tophat_flagstat");
	$pm->set_job_dependency($qid->{remove_duplicate});
	$pipeline->tophat->samtools_flagstat(
			sam          => "$pm->{dir}/$sample_id.mkdup.bam",
			output       => "$pm->{dir}/$sample_id.mkdup.flagstat",
			delete_input => 0,
		);

	
	####################################################################
	# more detailed QC
	####################################################################
	$pm->set_job_name("$sample_id"."_rnaseq_qc");
	$pm->set_job_dependency($qid->{remove_duplicate});
	$qid->{qc} = $pipeline->tophat->rnaseqqc(
		bam       => "$pm->{dir}/$sample_id.mkdup.bam",
		sample_id => $sample_id,
	);

	#####################################################################
	# counting
	####################################################################
	$pm->set_job_name("$sample_id"."_counting");
	$pm->set_job_dependency($qid->{remove_duplicate});
	$qid->{counting} = $pipeline->tophat->counting(
		bam => "$pm->{dir}/$sample_id.mkdup.bam",
		strand => $is_strand_specific,
		remove_duplicate => $remove_duplicate,
	);
}


1;
