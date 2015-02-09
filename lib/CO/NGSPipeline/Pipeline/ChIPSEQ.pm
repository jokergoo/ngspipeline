package CO::NGSPipeline::Pipeline::ChIPSEQ;

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
				 @_);
				 
	my $sample_id = $param{sample_id};
	my $r1 = $param{r1};
	my $r2 = $param{r2};
	my $library = $param{library};
	my $dir = $param{dir};
	
	my $qid = {};
	$qid->{flagstat} = [];
	$qid->{fastqc_r1} = [];
	$qid->{fastqc_r2} = [];
	my $sam_sort_list = [];
	
	for(my $i = 0; $i < scalar(@$r1); $i ++) {
		my $r1_fastq = $r1->[$i];
		my $r2_fastq = $r2->[$i];
		
	
		
		###################################################################
		# fastqc
		###################################################################
		$pm->set_job_name("$sample_id"."_chipseq_fastqc_r1_$i");
		$qid->{fastqc_r1}->[$i] = $pipeline->chipseq->fastqc(
			fastq      => $r1_fastq,
			output_dir => "$pm->{dir}/fastqc_r1_$i"
		);
		
		if($r2_fastq) {
			$pm->set_job_name("$sample_id"."_chipseq_fastqc_r2_$i");
			$qid->{fastqc_r2}->[$i] = $pipeline->chipseq->fastqc(
				fastq      => $r2_fastq,
				output_dir => "$pm->{dir}/fastqc_r2_$i"
			);
		}
	
		####################################################################
		# trim
		####################################################################
		$pm->set_job_name("$sample_id"."_chipseq_trimmed_$i");
		$qid->{trim} = $pipeline->chipseq->trim(
			fastq1  => $r1_fastq,
			fastq2  => $r2_fastq,
			output1 => "$pm->{dir}/$sample_id.r1.trimmed.$i.fastq.gz",
			output2 => "$pm->{dir}/$sample_id.r2.trimmed.$i.fastq.gz",
		);

				
		###################################################################
		# fastqc after trimming
		###################################################################
		$pm->set_job_name("$sample_id"."_chipseq_fastqc_r1_trimmed_$i");
		$pm->set_job_dependency($qid->{trim});
		$qid->{fastqc_r1_trimmed} = $pipeline->chipseq->fastqc(
			fastq      => "$pm->{dir}/$sample_id.r1.trimmed.$i.fastq.gz",
			output_dir => "$pm->{dir}/fastqc_r1_trimmed_$i",
		);

		if($r2_fastq) {
			$pm->set_job_name("$sample_id"."_chipseq_fastqc_r2_trimmed_$i");
			$pm->set_job_dependency($qid->{trim});
			$qid->{fastqc_r2_trimmed} = $pipeline->chipseq->fastqc(
				fastq      => "$pm->{dir}/$sample_id.r2.trimmed.$i.fastq.gz",
				output_dir => "$pm->{dir}/fastqc_r2_trimmed_$i",
			);
		}
	
		###################################################################
		# alignment
		###################################################################
		$pm->set_job_name("$sample_id"."_chipseq_alignment_r1_$i");
		$pm->set_job_dependency($qid->{fastqc_r1_trimmed});
		$qid->{alignment} = $pipeline->chipseq->bowtie_aln(
			fastq1        => "$pm->{dir}/$sample_id.r1.trimmed.$i.fastq.gz",
			fastq2        => "$pm->{dir}/$sample_id.r2.trimmed.$i.fastq.gz",
			output       => "$pm->{dir}/$sample_id.$i.sorted.bam",
			delete_input => 1,
		);
		
		
		####################################################################
		# flagstat
		####################################################################
		$pm->set_job_name("$sample_id"."_chipseq_flagstat_$i");
		$pm->set_job_dependency($qid->{alignment});
		$qid->{flagstat}->[$i] = $pipeline->chipseq->samtools_flagstat(
			sam          => "$pm->{dir}/$sample_id.$i.sorted.bam",
			output       => "$pm->{dir}/$sample_id.$i.flagstat",
			delete_input => 0
		);

		$sam_sort_list->[$i] = "$pm->{dir}/$sample_id.$i.sorted.bam";
	}

	########################################################################
	# merge, nodup
	########################################################################
	$pm->set_job_name("$sample_id"."_chipseq_merge_and_mkdup");
	$pm->set_job_dependency(@{$qid->{flagstat}});
	$qid->{remove_duplicate} = $pipeline->chipseq->merge_nodup(
		sam_list     => $sam_sort_list,
		output       => "$pm->{dir}/$sample_id.mkdup.bam",
		library      => $library,
		REMOVE_DUPLICATES => 'FALSE',
		delete_input => 1
	);
	

	####################################################################
	# flagstat
	####################################################################
	$pm->set_job_name("$sample_id"."_chipseq_flagstat_final");
	$pm->set_job_dependency($qid->{remove_duplicate});
	$pipeline->chipseq->samtools_flagstat(
		sam          => "$pm->{dir}/$sample_id.mkdup.bam",
		output       => "$pm->{dir}/$sample_id.mkdup.flagstat",
		delete_input => 0
	);	


	$pm->set_job_name("$sample_id"."_chipseq_data_transform");
	$pm->set_job_dependency($qid->{remove_duplicate});
	$qid->{data_transform} = $pipeline->chipseq->data_transform(
		sam          => "$pm->{dir}/$sample_id.mkdup.bam",
	);	
		
	$pm->set_job_name("$sample_id"."_chipseq_homer_qc");
	$pm->set_job_dependency($qid->{remove_duplicate});
	$qid->{homer_qc} = $pipeline->chipseq->homer_qc(
		sam          => "$pm->{dir}/$sample_id.mkdup.bam",
		output_dir   => $pm->{dir},
		sample_id    => $sample_id,
	);	
				
	if($sample_id ne "Pool_of_21_AK_samples") {
	
		$pm->set_job_name("$sample_id"."_chipseq_peak_calling");
		$pm->set_job_dependency($qid->{homer_qc});
		$qid->{peak_calling} = $pipeline->chipseq->peak_calling(
			sam          => "$pm->{dir}/$sample_id.mkdup.bam",
			output_dir   => $pm->{dir},
			sample_id    => $sample_id,
			homer_qc_dir => "$pm->{dir}/homer_qc",
		);	
	}
}

1;
