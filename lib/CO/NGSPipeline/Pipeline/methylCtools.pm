package CO::NGSPipeline::Pipeline::methylCtools;

use strict;
use CO::NGSPipeline::Tool::Config;
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
	my $sam_sort_list = [];
	for(my $i = 0; $i < scalar(@$r1); $i ++) {
		my $r1_fastq = $r1->[$i];
		my $r2_fastq = $r2->[$i];
			
		###################################################################
		# fastqc
		###################################################################
		$pm->set_job_name("$sample_id"."_methylctools_fastqc_r1_$i");
		$qid->{fastqc_r1} = $pipeline->methylctools->fastqc(
			fastq      => $r1_fastq,
			output_dir => "$pm->{dir}/fastqc_r1_$i"
		);

		$pm->set_job_name("$sample_id"."_methylctools_fastqc_r2_$i");
		$qid->{fastqc_r2} = $pipeline->methylctools->fastqc(
			fastq      => $r2_fastq,
			output_dir => "$pm->{dir}/fastqc_r2_$i"
		);

		####################################################################
		# trim
		####################################################################
		$pm->set_job_name("$sample_id"."_methylctools_trimmed_$i");
		$qid->{trim} = $pipeline->methylctools->trim(
			fastq1  => $r1_fastq,
			fastq2  => $r2_fastq,
			output1 => "$prefix1.trimmed.$i.fastq.gz",
			output2 => "$prefix2.trimmed.$i.fastq.gz",
		);	
				
		###################################################################
		# fastqc after trimming
		###################################################################
		$pm->set_job_name("$sample_id"."_methylctools_fastqc_r1_trimmed_$i");
		$pm->set_job_dependency($qid->{trim});
		$qid->{fastqc_r1_trimmed} = $pipeline->methylctools->fastqc(
			fastq        => "$prefix1.trimmed.$i.fastq.gz",
			output_dir   => "$pm->{dir}/fastqc_r1_trimmed_$i",
			delete_input => 0
		);

		$pm->set_job_name("$sample_id"."_methylctools_fastqc_r2_trimmed_$i");
		$pm->set_job_dependency($qid->{trim});
		$qid->{fastqc_r2_trimmed} = $pipeline->methylctools->fastqc(
			fastq => "$prefix2.trimmed.$i.fastq.gz",
			output_dir => "$pm->{dir}/fastqc_r2_trimmed_$i",
			delete_input => 0
		);
												 
		####################################################################
		# fqconv
		####################################################################
		$pm->set_job_name("$sample_id"."_methylctools_fqconv_ct_$i");
		$pm->set_job_dependency($qid->{trim});
		$qid->{fqconv_ct} = $pipeline->methylctools->fqconv(
			fastq        => "$prefix1.trimmed.$i.fastq.gz",
			output       => "$prefix1.trimmed.conv.$i.fastq.gz",
			which_pair   => 1,
			delete_input => 1
		);
				
		$pm->set_job_name("$sample_id"."_methylctools_fqconv_ga_$i");
		$pm->set_job_dependency($qid->{trim});
		$qid->{fqconv_ga} = $pipeline->methylctools->fqconv(
			fastq        => "$prefix2.trimmed.$i.fastq.gz",
			output       => "$prefix2.trimmed.conv.$i.fastq.gz",
			which_pair   => 2,
			delete_input => 1
		);
		
		####################################################################
		# align
		####################################################################
		my $use_convey = rand(1) > 0.8 ? 1 : 0;
		$pm->set_job_name("$sample_id"."_methylctools_alignment_ct_$i");
		$pm->set_job_dependency($qid->{fqconv_ct});
		$qid->{alignment_ct} = $pipeline->methylctools->bwa_aln(
			fastq        => "$prefix1.trimmed.conv.$i.fastq.gz",
			genome       => "$METHYLCTOOLS_GENOME_DIR/$METHYLCTOOLS_REF_GENOME_CONV",
			output       => "$prefix1.conv.$i.sai",
			delete_input => 0,
			use_convey   => $use_convey
		);
		
		$pm->set_job_name("$sample_id"."_methylctools_alignment_ga_$i");
		$pm->set_job_dependency($qid->{fqconv_ga});
		$qid->{alignment_ga} = $pipeline->methylctools->bwa_aln(
			fastq        => "$prefix2.trimmed.conv.$i.fastq.gz",
			genome       => "$METHYLCTOOLS_GENOME_DIR/$METHYLCTOOLS_REF_GENOME_CONV",
			output       => "$prefix2.conv.$i.sai",
			delete_input => 0,
			use_convey   => $use_convey
		);
				
		$pm->set_job_name("$sample_id"."_methylctools_sampe_$i");
		$pm->set_job_dependency($qid->{alignment_ct}, $qid->{alignment_ga});
		$qid->{sampe} = $pipeline->methylctools->sampe(
			aln1         => "$prefix1.conv.$i.sai",
			aln2         => "$prefix2.conv.$i.sai",
			fastq1       => "$prefix1.trimmed.conv.$i.fastq.gz",
			fastq2       => "$prefix2.trimmed.conv.$i.fastq.gz",
			genome       => "$METHYLCTOOLS_GENOME_DIR/$METHYLCTOOLS_REF_GENOME_CONV",
			output       => "$prefix1.conv.$i.bam",
			delete_input => 1
		);	
		
		####################################################################
		# bconv
		####################################################################
		$pm->set_job_name("$sample_id"."_methylctools_bconv_$i");
		$pm->set_job_dependency($qid->{sampe});
		$qid->{bconv} = $pipeline->methylctools->bconv(
			bam          => "$prefix1.conv.$i.bam",
			output       => "$prefix1.$i.bam",
			delete_input => 1
		);
		
		####################################################################
		# flagstat
		####################################################################
		$pm->set_job_name("$sample_id"."_methylctools_flagstat_$i");
		$pm->set_job_dependency($qid->{bconv});
		$qid->{flagstat} = $pipeline->methylctools->samtools_flagstat(
			sam          => "$prefix1.$i.bam",
			output       => "$prefix1.$i.flagstat",
			delete_input => 0
		);
		
		####################################################################
		# sort
		####################################################################
		$pm->set_job_name("$sample_id"."_methylctools_sort_$i");
		$pm->set_job_dependency($qid->{bconv});
		$qid->{sort_sam}->[$i] = $pipeline->methylctools->sort_sam(
			sam          => "$prefix1.$i.bam",
			output       => "$prefix1.sorted.$i.bam",
			delete_input => 1
		);
		$sam_sort_list->[$i] = "$prefix1.sorted.$i.bam";
	}
			
	########################################################################
	# merge, nodup
	########################################################################
	$pm->set_job_name("$sample_id"."_methylctools_merge_and_nodup");
	$pm->set_job_dependency(@{$qid->{sort_sam}});
	$qid->{remove_duplicate} = $pipeline->methylctools->merge_nodup(
		sam_list     => $sam_sort_list,
		output       => "$prefix1.nodup.bam", 
		library      => $library,
		delete_input => 1
	);

	########################################################################
	# insert size
	########################################################################
	$pm->set_job_name("$sample_id"."_methylctools_insertsize");
	$pm->set_job_dependency($qid->{remove_duplicate});
	$qid->{insertsize} = $pipeline->methylctools->picard_insertsize(
		sam => "$prefix1.nodup.bam"
	);

	########################################################################
	# lambda conversion
	########################################################################
	$pm->set_job_name("$sample_id"."_methylctools_lambda_conversion");
	$pm->set_job_dependency($qid->{remove_duplicate});
	$qid->{lambda_conversion} = $pipeline->methylctools->lambda_conversion(
		bam    => "$prefix1.nodup.bam",
		output => "$prefix1.lambda.conversion.txt"
	);
													 
	########################################################################
	# methylation calling
	########################################################################
	if($no_bissnp) {
		$pm->set_job_name("$sample_id"."_methylctools_bcall");
		$pm->set_job_dependency($qid->{remove_duplicate});
		$qid->{bcall} = $pipeline->methylctools->bcall(
			bam    => "$prefix1.nodup.bam",
			output => "$prefix1.call.gz"
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
												
		$pm->set_job_name("$sample_id"."_methylctools_QC");
		$pm->set_job_dependency($qid->{methy_calling});
		$qid->{qc} = $pipeline->methylctools->bsqc(
			dir => $pm->{dir},
			tool => "methylctools",
			sample => "$sample_id",
		);
		
		for my $chr (map {"chr$_"} (1..22, "X", "Y")) {
			$pm->set_job_name("$sample_id"."_bsseq_RData_$chr");
			$pm->set_job_dependency($qid->{methy_calling});
			$qid->{RData} = $pipeline->bsmooth->RData(
				input   => "$prefix1.nodup.cpg.filtered.CG.bedgraph",
				sample_id  => $sample_id,
				chr        => $chr,
				output     => "$pm->{dir}/bsseq.$chr.RData",
			);
		}
	}
}

1;
