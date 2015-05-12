package CO::NGSPipeline::Pipeline::STAR;

use strict;
use File::Basename;
use base qw/CO::NGSPipeline 
            CO::NGSPipeline::Tool::RNAseq::STAR/;

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
		$pm->set_job_name("$sample_id"."_star_fastqc_r1_$i");
		$qid->{fastqc_r1} = $pipeline->star->fastqc(
				fastq      => $r1_fastq,
				output_dir => "$pm->{dir}/fastqc_r1_$i"
		);

		if($r2_fastq) {
			$pm->set_job_name("$sample_id"."_star_fastqc_r2_$i");
			$qid->{fastqc_r2} = $pipeline->star->fastqc(
					fastq      => $r2_fastq,
					output_dir => "$pm->{dir}/fastqc_r2_$i"
			);
		}

		##################################################################
		# alignment
		##################################################################
		$pm->set_job_name("$sample_id"."_star_align_$i");
		$pm->set_job_dependency($qid->{trim});
		$qid->{alignment} = $pipeline->star->align(
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
		$pm->set_job_name("$sample_id"."_star_flagstat_$i");
		$pm->set_job_dependency($qid->{alignment});
		$qid->{flagstat}->[$i] = $pipeline->star->samtools_flagstat(
			sam          => "$pm->{dir}/$sample_id.$i.bam",
			output       => "$pm->{dir}/$sample_id.$i.flagstat",
			delete_input => 0,
		);
		
		$sam_sort_list->[$i] = "$pm->{dir}/$sample_id.$i.bam";
	}
	
	$pm->set_job_name("$sample_id"."_star_merge_and_mkdup");
	$pm->set_job_dependency(@{$qid->{flagstat}});
	$qid->{remove_duplicate} = $pipeline->star->merge_nodup(
			sam_list     => $sam_sort_list,
			output       => "$pm->{dir}/$sample_id.mkdup.bam",
			delete_input => 1,
			REMOVE_DUPLICATES => 'FALSE',
	);

	$pm->set_job_name("$sample_id"."_star_flagstat");
	$pm->set_job_dependency($qid->{remove_duplicate});
	$pipeline->star->samtools_flagstat(
			sam          => "$pm->{dir}/$sample_id.mkdup.bam",
			output       => "$pm->{dir}/$sample_id.mkdup.flagstat",
			delete_input => 0,
		);

	
	####################################################################
	# more detailed QC
	####################################################################
	$pm->set_job_name("$sample_id"."_rnaseq_qc");
	$pm->set_job_dependency($qid->{remove_duplicate});
	$qid->{qc} = $pipeline->star->rnaseqqc(
		bam       => "$pm->{dir}/$sample_id.mkdup.bam",
		sample_id => $sample_id,
	);

	####################################################################
	# sort by name
	####################################################################
	$pm->set_job_name("$sample_id"."_namesort");
	$pm->set_job_dependency($qid->{remove_duplicate});
	$qid->{namesort} = $pipeline->star->sort_sam(
		sam     => "$pm->{dir}/$sample_id.mkdup.bam",
		output  => "$pm->{dir}/$sample_id.mkdup.namesorted.bam",
		sort_by => "queryname",
		delete_input => 0,
	);

	#####################################################################
	# counting
	####################################################################
	$pm->set_job_name("$sample_id"."_counting");
	$pm->set_job_dependency($qid->{qc}, $qid->{namesort});
	$qid->{counting} = $pipeline->star->counting(
		bam => "$pm->{dir}/$sample_id.mkdup.namesorted.bam",
		sample_id => $sample_id,
	);
}

1;
