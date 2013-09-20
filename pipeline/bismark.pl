

# importing from mother script:
## for every sample
#     $pipeline : CO::NGSPipeline object
#     $r1       : pair 1, array reference
#     $r2       : pair 2, array reference
#     $prefix1  : prefix for pair 1 to output
#     $prefix2  : prefix for pari 2 to output
#     $library  : library class for lanes of the sample

my $qid = {};
$qid->{sort_sam} = [];
my $sam_sorted_list = [];
for(my $i = 0; $i < scalar(@$r1); $i ++) {

	my $r1_fastq = $r1->[$i];
	my $r2_fastq = $r2->[$i];
	
	###################################################################
	# fastqc
	###################################################################
	$pipeline->set_job_name("$sample_id"."_bismark_fastqc_r1_$i");
	$qid->{fastqc_r1} = $pipeline->bismark->fastqc(
		fastq      => $r1_fastq,
		output_dir => "$pipeline->{dir}/fastqc_r1_$i"
	);

	$pipeline->set_job_name("$sample_id"."_bismark_fastqc_r2_$i");
	$qid->{fastqc_r2} = $pipeline->bismark->fastqc(
		fastq      => $r2_fastq,
		output_dir => "$pipeline->{dir}/fastqc_r2_$i"
	);
			
	####################################################################
	# trim
	####################################################################
	$pipeline->set_job_name("$sample_id"."_bismark_trimmed_$i");
	$qid->{trim} = $pipeline->bismark->trim(
		fastq1  => $r1_fastq,
		fastq2  => $r2_fastq,
		output1 => "$prefix1.trimmed.$i.fastq.gz",
		output2 => "$prefix2.trimmed.$i.fastq.gz"
	);
	
	###################################################################
	# fastqc after trimming
	###################################################################
	$pipeline->set_job_name("$sample_id"."_bismark_fastqc_r1_trimmed_$i");
	$pipeline->set_job_dependency($qid->{trim});
	$qid->{fastqc_r1_trimmed} = $pipeline->bismark->fastqc(
		fastq      => "$prefix1.trimmed.$i.fastq.gz",
		output_dir => "$pipeline->{dir}/fastqc_r1_trimmed_$i"
	);

	$pipeline->set_job_name("$sample_id"."_bismark_fastqc_r2_trimmed_$i");
	$pipeline->set_job_dependency($qid->{trim});
	$qid->{fastqc_r2_trimmed} = $pipeline->bismark->fastqc(
		fastq      => "$prefix2.trimmed.$i.fastq.gz",
		output_dir => "$pipeline->{dir}/fastqc_r2_trimmed_$i"
	);
	 
	####################################################################
	# align
	####################################################################
	$pipeline->set_job_name("$sample_id"."_bismark_alignment_$i");
	$pipeline->set_job_dependency($qid->{trim});
	$qid->{alignment} = $pipeline->bismark->align(
		fastq1       => "$prefix1.trimmed.$i.fastq.gz",
		fastq2       => "$prefix2.trimmed.$i.fastq.gz",
		output       => "$prefix1.$i.bam",  # must be .bam
		delete_input => 1
	);
	
	####################################################################
	# flagstat
	####################################################################
	$pipeline->set_job_name("$sample_id"."_bismark_flagstat_$i");
	$pipeline->set_job_dependency($qid->{alignment});
	$qid->{flagstat} = $pipeline->bismark->samtools_flagstat(
		sam          => "$prefix1.$i.bam",
		output       => "$prefix1.$i.flagstat",
		delete_input => 0
	);

	####################################################################
	# sort
	####################################################################
	$pipeline->set_job_name("$sample_id"."_bismark_sort_sam_$i");
	$pipeline->set_job_dependency($qid->{alignment});
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
$pipeline->set_job_name("$sample_id"."_bismark_merge_and_nodup");
$pipeline->set_job_dependency(@{$qid->{sort_sam}});
$qid->{remove_duplicate} = $pipeline->bismark->merge_nodup(
	sam_list     => $sam_sorted_list,
	output       => "$prefix1.nodup.bam", 
	library      => $library,
	delete_input => 1
);
										   
########################################################################
# insert size
########################################################################
$pipeline->set_job_name("$sample_id"."_bismark_insertsize");
$pipeline->set_job_dependency($qid->{remove_duplicate});
$qid->{insertsize} = $pipeline->bismark->picard_insertsize(
	sam => "$prefix1.nodup.bam"
);
		
########################################################################
# lambda conversion
########################################################################
$pipeline->set_job_name("$sample_id"."_bismark_lambda_conversion");
$pipeline->set_job_dependency($qid->{remove_duplicate});
$qid->{lambda_conversion} = $pipeline->bismark->lambda_conversion(
	bam    => "$prefix1.nodup.bam",
	output => "$prefix1.lambda.conversion.txt"
);
		
if($no_bissnp) {
	####################################################################
	# sort by queryname
	####################################################################
	$pipeline->set_job_name("$sample_id"."_bismark_sort_sam_by_queryname");
	$pipeline->set_job_dependency($qid->{remove_duplicate});
	$qid->{sort_sam_query_name} = $pipeline->bismark->sort_sam(
		sam          => "$prefix1.nodup.bam",
		output       => "$prefix1.nodup.queryname_sorted.bam",
		sort_by      => "queryname",
		delete_input => 0
	);

	########################################################################
	# methylation calling
	########################################################################
	$pipeline->set_job_name("$sample_id"."_bismark_methylation_calling");
	$pipeline->set_job_dependency($qid->{sort_sam_query_name});
	$qid->{call_methylation} = $pipeline->bismark->call_methylation(
		sam          => "$prefix1.nodup.queryname_sorted.bam",
		delete_input => 1
	);
	
} else {
	########################################################################
	# bissnp methylation calling
	########################################################################
	$pipeline->set_job_name("$sample_id"."_bissnp_methylation_calling");
	$pipeline->set_job_dependency($qid->{remove_duplicate});
	$qid->{methy_calling} = $pipeline->bissnp->call_methylation(
		bam => "$prefix1.nodup.bam"
	);
											
	$pipeline->set_job_name("$sample_id"."_bismark_QC");
	$pipeline->set_job_dependency($qid->{methy_calling});
	$qid->{qc} = $pipeline->bismark->bsqc(
		dir    => $pipeline->{dir},
		tool   => "bismark",
		sample => "$sample_id",
		base_dir => $SCRIPT_DIR,
	);
	
	for my $chr (map {"chr$_"} (1..22, "X", "Y")) {
		$pipeline->set_job_name("$sample_id"."_bsseq_RData_$chr");
		$pipeline->set_job_dependency($qid->{methy_calling});
		$qid->{RData} = $pipeline->bismark->RData(
			bedgraph   => "$prefix1.nodup.cpg.filtered.CG.bedgraph",
			sample_id  => $sample_id,
			chr        => $chr,
		);
	}
}


1;
