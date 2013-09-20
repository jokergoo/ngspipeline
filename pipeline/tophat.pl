# importing from mother script:
## for every sample
#     $pipeline : CO::NGSPipeline object
#     $r1       : pair 1, array reference
#     $r2       : pair 2, array reference
#     $prefix1  : prefix for pair 1 to output
#     $prefix2  : prefix for pari 2 to output
#     $library  : library class for lanes of the sample

my $qid;

my $r1_fastq = $r1->[0];
my $r2_fastq = $r2->[0];

###################################################################
# fastqc
###################################################################
$pipeline->set_job_name("$sample_id"."_tophat_fastqc_r1_$i");
$qid->{fastqc_r1} = $pipeline->tophat->fastqc(
		fastq      => $r1_fastq,
		output_dir => "$pipeline->{dir}/fastqc_r1_$i"
);

$pipeline->set_job_name("$sample_id"."_tophat_fastqc_r2_$i");
$qid->{fastqc_r2} = $pipeline->tophat->fastqc(
		fastq      => $r2_fastq,
		output_dir => "$pipeline->{dir}/fastqc_r2_$i"
);

####################################################################
# trim
####################################################################
$pipeline->set_job_name("$sample_id"."_tophat_trimmed_$i");
$qid->{trim} = $pipeline->tophat->trim(
	fastq1  => $r1_fastq,
	fastq2  => $r2_fastq,
	output1 => "$prefix1.trimmed.$i.fastq.gz",
	output2 => "$prefix2.trimmed.$i.fastq.gz",
	polya   => 1,
);	
			
###################################################################
# fastqc after trimming
###################################################################
$pipeline->set_job_name("$sample_id"."_tophat_fastqc_r1_trimmed_$i");
$pipeline->set_job_dependency($qid->{trim});
$qid->{fastqc_r1_trimmed} = $pipeline->tophat->fastqc(
	fastq        => "$prefix1.trimmed.$i.fastq.gz",
	output_dir   => "$pipeline->{dir}/fastqc_r1_trimmed_$i",
	delete_input => 0
);

$pipeline->set_job_name("$sample_id"."_tophat_fastqc_r2_trimmed_$i");
$pipeline->set_job_dependency($qid->{trim});
$qid->{fastqc_r2_trimmed} = $pipeline->tophat->fastqc(
	fastq        => "$prefix2.trimmed.$i.fastq.gz",
	output_dir   => "$pipeline->{dir}/fastqc_r2_trimmed_$i",
	delete_input => 0
);

##################################################################
# alignment
##################################################################
$pipeline->set_job_name("$sample_id"."_tophat_align");
$pipeline->set_job_dependency($qid->{trim});
$qid->{align} = $pipeline->tophat->align(
	fastq1 => "$prefix1.trimmed.$i.fastq.gz",
	fastq2 => "$prefix2.trimmed.$i.fastq.gz",
	output => "$pipeline->{dir}/$sample_id.bam",
	sample_id => $sample_id,
	strand => $is_strand_specific,
	delete_input => 1
);

####################################################################
# more detailed QC
####################################################################
$pipeline->set_job_name("$sample_id"."_rnaseq_qc");
$pipeline->set_job_dependency($qid->{align});
$qid->{qc} = $pipeline->tophat->rnaseqqc(
	bam       => "$pipeline->{dir}/$sample_id.bam",
	sample_id => $sample_id
);

#####################################################################
# counting and RPKM
####################################################################
$pipeline->set_job_name("$sample_id"."_counting");
$pipeline->set_job_dependency($qid->{align});
$qid->{counting} = $pipeline->tophat->counting(
	bam => "$pipeline->{dir}/$sample_id.bam",
	strand => $is_strand_specific
);

1;
