package CO::NGSPipeline::Pipeline::BSMAP;

use strict;
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
	$qid->{alignment} = [];
	my $sam_sort_list = [];
	for(my $i = 0; $i < scalar(@$r1); $i ++) {
		my $r1_fastq = $r1->[$i];
		my $r2_fastq = $r2->[$i];
		
		
		###################################################################
		# fastqc
		###################################################################
		$pm->set_job_name("$sample_id"."_bsmap_fastqc_r1_$i");
		$qid->{fastqc_r1} = $pipeline->bsmap->fastqc(
			fastq      => $r1_fastq,
			output_dir => "$pm->{dir}/fastqc_r1_$i"
		);

		$pm->set_job_name("$sample_id"."_bsmap_fastqc_r2_$i");
		$qid->{fastqc_r2} = $pipeline->bsmap->fastqc(
			fastq      => $r2_fastq,
			output_dir => "$pm->{dir}/fastqc_r2_$i"
		);
				
		####################################################################
		# trim
		####################################################################
		$pm->set_job_name("$sample_id"."_bsmap_trimmed_$i");
		$qid->{trim} = $pipeline->bsmap->trim(
			fastq1  => $r1_fastq,
			fastq2  => $r2_fastq,
			output1 => "$prefix1.trimmed.$i.fastq.gz",
			output2 => "$prefix2.trimmed.$i.fastq.gz",
		);
				
		###################################################################
		# fastqc after trimming
		###################################################################
		$pm->set_job_name("$sample_id"."_bsmap_fastqc_r1_trimmed_$i");
		$pm->set_job_dependency($qid->{trim});
		$qid->{fastqc_r1_trimmed} = $pipeline->bsmap->fastqc(
			fastq      => "$prefix1.trimmed.$i.fastq.gz",
			output_dir => "$pm->{dir}/fastqc_r1_trimmed_$i"
		);

		$pm->set_job_name("$sample_id"."_bsmap_fastqc_r2_trimmed_$i");
		$pm->set_job_dependency($qid->{trim});
		$qid->{fastqc_r2_trimmed} = $pipeline->bsmap->fastqc(
			fastq      => "$prefix2.trimmed.$i.fastq.gz",
			output_dir => "$pm->{dir}/fastqc_r2_trimmed_$i"
		);
												 
		###################################################################
		# alignment
		###################################################################
		$pm->set_job_name("$sample_id"."_bsmap_alignment_$i");
		$pm->set_job_dependency($qid->{trim});
		$qid->{alignment}->[$i] = $pipeline->bsmap->align(
			fastq1       => "$prefix1.trimmed.$i.fastq.gz",
			fastq2       => "$prefix2.trimmed.$i.fastq.gz",
			output       => "$prefix1.$i.sorted.bam",
			delete_input => 1,
		);
				
		####################################################################
		# flagstat
		####################################################################
		$pm->set_job_name("$sample_id"."_bsmap_flagstat_$i");
		$pm->set_job_dependency($qid->{alignment}->[$i]);
		$qid->{flagstat} = $pipeline->bsmap->samtools_flagstat(
			sam          => "$prefix1.$i.sorted.bam",
			output       => "$prefix1.$i.flagstat",
			delete_input => 0
		);

		$sam_sort_list->[$i] = "$prefix1.$i.sorted.bam";
	}

	########################################################################
	# merge, nodup
	########################################################################
	$pm->set_job_name("$sample_id"."_bsmap_merge_and_nodup");
	$pm->set_job_dependency(@{$qid->{alignment}});
	$qid->{remove_duplicate} = $pipeline->bsmap->merge_nodup(
		sam_list     => $sam_sort_list,
		output       => "$prefix1.nodup.bam",
		library      => $library,
		delete_input => 1
	);
									   
	########################################################################
	# insert size
	########################################################################
	$pm->set_job_name("$sample_id"."_bsmap_insertsize");
	$pm->set_job_dependency($qid->{remove_duplicate});
	$qid->{insertsize} = $pipeline->bsmap->picard_insertsize(
		sam => "$prefix1.nodup.bam"
	);

	########################################################################
	# lambda conversion
	########################################################################
	$pm->set_job_name("$sample_id"."_bsmap_lambda_conversion");
	$pm->set_job_dependency($qid->{remove_duplicate});
	$qid->{lambda_conversion} = $pipeline->bsmap->lambda_conversion(
		bam    => "$prefix1.nodup.bam",
		output => "$prefix1.lambda.conversion.txt"
	);
				
	########################################################################
	# call methylation, default methylation calling by BSMAP
	########################################################################
	if($no_bissnp) {
		$pm->set_job_name("$sample_id"."_bsmap_methylation_calling");
		$pm->set_job_dependency($qid->{remove_duplicate});
		$qid->{methy_calling} = $pipeline->bsmap->call_methylation(
			sam    => "$prefix1.nodup.bam",
			output => "$prefix1.methylcall.txt"
		);
	} else {
	########################################################################
	# bissnp methylation calling
	########################################################################
		$pm->set_job_name("$sample_id"."_bissnp_methylation_calling");
		$pm->set_job_dependency($qid->{remove_duplicate});
		$qid->{methy_calling} = $pipeline->bissnp->call_methylation(
			bam => "$prefix1.nodup.bam",
		);
												
		$pm->set_job_name("$sample_id"."_bsmap_QC");
		$pm->set_job_dependency($qid->{methy_calling});
		$qid->{qc} = $pipeline->bsmap->bsqc(
			dir    => $pm->{dir},
			tool   => "bsmap",
			sample => "$sample_id",
			base_dir => $SCRIPT_DIR,
		);
		
		for my $chr (map {"chr$_"} (1..22, "X", "Y")) {
			$pm->set_job_name("$sample_id"."_bsseq_RData_$chr");
			$pm->set_job_dependency($qid->{methy_calling});
			$qid->{RData} = $pipeline->bsmap->RData(
				bedgraph   => "$prefix1.nodup.cpg.filtered.CG.bedgraph",
				sample_id  => $sample_id,
				chr        => $chr,
			);
		}
	}
}

1;
