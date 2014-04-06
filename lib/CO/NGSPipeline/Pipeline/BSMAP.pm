package CO::NGSPipeline::Pipeline::BSMAP;

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
	             r1        => undef,
				 r2        => undef,
				 library   => undef,
				 dir       => undef,
				 no_bissnp => 0,
				 species   => 'human',
				 @_);
				 
	my $sample_id = $param{sample_id};
	my $r1        = $param{r1};
	my $r2        = $param{r2};
	my $library   = $param{library};
	my $dir       = $param{dir};
	my $no_bissnp = $param{no_bissnp};
	my $species   = $param{species};

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
		$pm->set_job_name("$sample_id"."_bsmap_fastqc_r1_$i");
		$qid->{fastqc_r1}->[$i] = $pipeline->bsmap->fastqc(
			fastq      => $r1_fastq,
			output_dir => "$pm->{dir}/fastqc_r1_$i"
		);

		$pm->set_job_name("$sample_id"."_bsmap_fastqc_r2_$i");
		$qid->{fastqc_r2}->[$i] = $pipeline->bsmap->fastqc(
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
			output1 => "$pm->{dir}/$sample_id.r1.trimmed.$i.fastq.gz",
			output2 => "$pm->{dir}/$sample_id.r2.trimmed.$i.fastq.gz",
		);
				
		###################################################################
		# fastqc after trimming
		###################################################################
		$pm->set_job_name("$sample_id"."_bsmap_fastqc_r1_trimmed_$i");
		$pm->set_job_dependency($qid->{trim});
		$qid->{fastqc_r1_trimmed} = $pipeline->bsmap->fastqc(
			fastq      => "$pm->{dir}/$sample_id.r1.trimmed.$i.fastq.gz",
			output_dir => "$pm->{dir}/fastqc_r1_trimmed_$i"
		);

		$pm->set_job_name("$sample_id"."_bsmap_fastqc_r2_trimmed_$i");
		$pm->set_job_dependency($qid->{trim});
		$qid->{fastqc_r2_trimmed} = $pipeline->bsmap->fastqc(
			fastq      => "$pm->{dir}/$sample_id.r2.trimmed.$i.fastq.gz",
			output_dir => "$pm->{dir}/fastqc_r2_trimmed_$i"
		);
												 
		###################################################################
		# alignment
		###################################################################
		$pm->set_job_name("$sample_id"."_bsmap_alignment_$i");
		$pm->set_job_dependency($qid->{fastqc_r1_trimmed}, $qid->{fastqc_r2_trimmed});
		$qid->{alignment} = $pipeline->bsmap->align(
			fastq1       => "$pm->{dir}/$sample_id.r1.trimmed.$i.fastq.gz",
			fastq2       => "$pm->{dir}/$sample_id.r2.trimmed.$i.fastq.gz",
			output       => "$pm->{dir}/$sample_id.$i.sorted.bam",
			delete_input => 1,
			species      => $species,
		);
				
		####################################################################
		# flagstat
		####################################################################
		$pm->set_job_name("$sample_id"."_bsmap_flagstat_$i");
		$pm->set_job_dependency($qid->{alignment});
		$qid->{flagstat}->[$i] = $pipeline->bsmap->samtools_flagstat(
			sam          => "$pm->{dir}/$sample_id.$i.sorted.bam",
			output       => "$pm->{dir}/$sample_id.$i.flagstat",
			delete_input => 0,
		);

		$sam_sort_list->[$i] = "$pm->{dir}/$sample_id.$i.sorted.bam";
	}

	########################################################################
	# merge, nodup
	########################################################################
	$pm->set_job_name("$sample_id"."_bsmap_merge_and_nodup");
	$pm->set_job_dependency(@{$qid->{flagstat}});
	$qid->{remove_duplicate} = $pipeline->bsmap->merge_nodup(
		sam_list     => $sam_sort_list,
		output       => "$pm->{dir}/$sample_id.nodup.bam",
		library      => $library,
		delete_input => 1,
	);
									   
	########################################################################
	# insert size
	########################################################################
	$pm->set_job_name("$sample_id"."_bsmap_insertsize");
	$pm->set_job_dependency($qid->{remove_duplicate});
	$qid->{insertsize} = $pipeline->bsmap->picard_insertsize(
		sam => "$pm->{dir}/$sample_id.nodup.bam",
	);

	########################################################################
	# lambda conversion
	########################################################################
	$pm->set_job_name("$sample_id"."_bsmap_lambda_conversion");
	$pm->set_job_dependency($qid->{remove_duplicate});
	$qid->{lambda_conversion} = $pipeline->bsmap->lambda_conversion(
		bam    => "$pm->{dir}/$sample_id.nodup.bam",
		output => "$pm->{dir}/$sample_id.lambda.conversion.txt",
		species => $species,
	);
				
	########################################################################
	# call methylation, default methylation calling by BSMAP
	########################################################################
	if($no_bissnp) {
		$pm->set_job_name("$sample_id"."_bsmap_methylation_calling");
		$pm->set_job_dependency($qid->{insertsize}, $qid->{lambda_conversion});
		$qid->{methy_calling} = $pipeline->bsmap->call_methylation(
			sam    => "$pm->{dir}/$sample_id.nodup.bam",
			output => "$pm->{dir}/$sample_id.methylcall.txt",
			species => $species,
		);
	} else {
	########################################################################
	# bissnp methylation calling
	########################################################################
		$pm->set_job_name("$sample_id"."_bissnp_methylation_calling");
		$pm->set_job_dependency($qid->{insertsize}, $qid->{lambda_conversion});
		$qid->{methy_calling} = $pipeline->bissnp->call_methylation(
			bam => "$pm->{dir}/$sample_id.nodup.bam",
			species => $species,
		);
												
		$pm->set_job_name("$sample_id"."_bsmap_QC");
		$pm->set_job_dependency(@{$qid->{fastqc_r1}}, @{$qid->{fastqc_r2}}, $qid->{methy_calling});
		$qid->{qc} = $pipeline->bsmap->bsqc(
			dir    => $pm->{dir},
			tool   => "bsmap",
			sample => "$sample_id",
		);
		
		for my $chr (map {"chr$_"} (1..22, "X", "Y")) {
			$pm->set_job_name("$sample_id"."_bsseq_RData_$chr");
			$pm->set_job_dependency($qid->{methy_calling});
			$qid->{RData} = $pipeline->bsmooth->RData(
				input   => "$pm->{dir}/$sample_id.nodup.cpg.filtered.CG.bedgraph",
				sample_id  => $sample_id,
				chr        => $chr,
				output     => "$pm->{dir}/bsseq.$chr.RData",
			);
		}
	}
}

1;
