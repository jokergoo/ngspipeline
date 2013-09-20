
my $qid;

my $r1_fastq = $r1->[0];
my $r2_fastq = $r2->[0];
			
$pipeline->set_job_name("$sample_id"."_defuse");
$qid = $pipeline->genefusion->defuse(
	fastq1 => $r1_fastq,
	fastq2 => $r2_fastq
);

1;
