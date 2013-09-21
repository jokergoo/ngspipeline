package CO::NGSPipeline::Pipeline::Bismark;

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
				 no_bissnp => 0,
				 @_);
				 
	my $sample_id = $param{sample_id};
	my $r1 = $param{r1};
	my $r2 = $param{r2};
	my $library = $param{library};
	my $dir = $param{dir};
	my $no_bissnp = $param{no_bissnp};
	
	# prefix means absolute path without fast/fq or fast.gz/fq.gz
	my $prefix1 = basename($r1->[0]);
	$prefix1 =~s/\.(fq|fastq)(\.gz)?$//;
	$prefix1 = "$pm->{dir}/$prefix1";
	my $prefix2 = basename($r2->[0]);
	$prefix2 =~s/\.(fq|fastq)(\.gz)?$//;
	$prefix2 = "$pm->{dir}/$prefix2";
	
	
	my $qid = {};
	$qid->{sort_sam} = [];
	my $sam_sorted_list = [];
	for(my $i = 0; $i < scalar(@$r1); $i ++) {

		my $r1_fastq = $r1->[$i];
		my $r2_fastq = $r2->[$i];
		
		###################################################################
		# fastqc
		###################################################################
		$pm->set_job_name("$sample_id"."_bismark_fastqc_r1_$i");
		$qid->{fastqc_r1} = $pipeline->bismark->fastqc(
			fastq      => $r1_fastq,
			output_dir => "$pm->{dir}/fastqc_r1_$i"
		);

		$pm->set_job_name("$sample_id"."_bismark_fastqc_r2_$i");
		$qid->{fastqc_r2} = $pipeline->bismark->fastqc(
			fastq      => $r2_fastq,
			output_dir => "$pm->{dir}/fastqc_r2_$i"
		);
				
		####################################################################
		# trim
		####################################################################
		$pm->set_job_name("$sample_id"."_bismark_trimmed_$i");
		$qid->{trim} = $pipeline->bismark->trim(
			fastq1  => $r1_fastq,
			fastq2  => $r2_fastq,
			output1 => "$prefix1.trimmed.$i.fastq.gz",
			output2 => "$prefix2.trimmed.$i.fastq.gz"
		);
		
		###################################################################
		# fastqc after trimming
		###################################################################
		$pm->set_job_name("$sample_id"."_bismark_fastqc_r1_trimmed_$i");
		$pm->set_job_dependency($qid->{trim});
		$qid->{fastqc_r1_trimmed} = $pipeline->bismark->fastqc(
			fastq      => "$prefix1.trimmed.$i.fastq.gz",
			output_dir => "$pm->{dir}/fastqc_r1_trimmed_$i"
		);

		$pm->set_job_name("$sample_id"."_bismark_fastqc_r2_trimmed_$i");
		$pm->set_job_dependency($qid->{trim});
		$qid->{fastqc_r2_trimmed} = $pipeline->bismark->fastqc(
			fastq      => "$prefix2.trimmed.$i.fastq.gz",
			output_dir => "$pm->{dir}/fastqc_r2_trimmed_$i"
		);
		 
		####################################################################
		# align
		####################################################################
		$pm->set_job_name("$sample_id"."_bismark_alignment_$i");
		$pm->set_job_dependency($qid->{trim});
		$qid->{alignment} = $pipeline->bismark->align(
			fastq1       => "$prefix1.trimmed.$i.fastq.gz",
			fastq2       => "$prefix2.trimmed.$i.fastq.gz",
			output       => "$prefix1.$i.bam",  # must be .bam
			delete_input => 1
		);
		
		####################################################################
		# flagstat
		####################################################################
		$pm->set_job_name("$sample_id"."_bismark_flagstat_$i");
		$pm->set_job_dependency($qid->{alignment});
		$qid->{flagstat} = $pipeline->bismark->samtools_flagstat(
			sam          => "$prefix1.$i.bam",
			output       => "$prefix1.$i.flagstat",
			delete_input => 0
		);

		####################################################################
		# sort
		####################################################################
		$pm->set_job_name("$sample_id"."_bismark_sort_sam_$i");
		$pm->set_job_dependency($qid->{alignment});
		$qid->{sort_sam}->[$i] = $pipeline->bismark->sort_sam(
			sam          => "$prefix1.$i.bam",
			output       => "$prefix1.$i.sorted.bam", 
			delete_input => 1
		);
		$sam_sorted_list->[$i] = "$prefix1.$i.sorted.bam";
		
	}
			
	########################################################################
	# merge, nodup
	########################################################################
	$pm->set_job_name("$sample_id"."_bismark_merge_and_nodup");
	$pm->set_job_dependency(@{$qid->{sort_sam}});
	$qid->{remove_duplicate} = $pipeline->bismark->merge_nodup(
		sam_list     => $sam_sorted_list,
		output       => "$prefix1.nodup.bam", 
		library      => $library,
		delete_input => 1
	);
											   
	########################################################################
	# insert size
	########################################################################
	$pm->set_job_name("$sample_id"."_bismark_insertsize");
	$pm->set_job_dependency($qid->{remove_duplicate});
	$qid->{insertsize} = $pipeline->bismark->picard_insertsize(
		sam => "$prefix1.nodup.bam"
	);
			
	########################################################################
	# lambda conversion
	########################################################################
	$pm->set_job_name("$sample_id"."_bismark_lambda_conversion");
	$pm->set_job_dependency($qid->{remove_duplicate});
	$qid->{lambda_conversion} = $pipeline->bismark->lambda_conversion(
		bam    => "$prefix1.nodup.bam",
		output => "$prefix1.lambda.conversion.txt"
	);
			
	if($no_bissnp) {
		####################################################################
		# sort by queryname
		####################################################################
		$pm->set_job_name("$sample_id"."_bismark_sort_sam_by_queryname");
		$pm->set_job_dependency($qid->{remove_duplicate});
		$qid->{sort_sam_query_name} = $pipeline->bismark->sort_sam(
			sam          => "$prefix1.nodup.bam",
			output       => "$prefix1.nodup.queryname_sorted.bam",
			sort_by      => "queryname",
			delete_input => 0
		);

		########################################################################
		# methylation calling
		########################################################################
		$pm->set_job_name("$sample_id"."_bismark_methylation_calling");
		$pm->set_job_dependency($qid->{sort_sam_query_name});
		$qid->{call_methylation} = $pipeline->bismark->call_methylation(
			sam          => "$prefix1.nodup.queryname_sorted.bam",
			delete_input => 1
		);
		
	} else {
		########################################################################
		# bissnp methylation calling
		########################################################################
		$pm->set_job_name("$sample_id"."_bissnp_methylation_calling");
		$pm->set_job_dependency($qid->{remove_duplicate});
		$qid->{methy_calling} = $pipeline->bissnp->call_methylation(
			bam => "$prefix1.nodup.bam"
		);
												
		$pm->set_job_name("$sample_id"."_bismark_QC");
		$pm->set_job_dependency($qid->{methy_calling});
		$qid->{qc} = $pipeline->bismark->bsqc(
			dir    => $pipeline->{dir},
			tool   => "bismark",
			sample => "$sample_id",
		);
		
		for my $chr (map {"chr$_"} (1..22, "X", "Y")) {
			$pm->set_job_name("$sample_id"."_bsseq_RData_$chr");
			$pm->set_job_dependency($qid->{methy_calling});
			$qid->{RData} = $pipeline->bismark->RData(
				bedgraph   => "$prefix1.nodup.cpg.filtered.CG.bedgraph",
				sample_id  => $sample_id,
				chr        => $chr,
			);
		}
	}

}

1;
