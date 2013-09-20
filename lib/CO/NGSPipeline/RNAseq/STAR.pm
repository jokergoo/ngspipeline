package CO::NGSPipeline::RNAseq::STAR;

use strict;
use CO::NGSPipeline::RNAseq::Config;
use CO::NGSPipeline::Utils;

use base qw(CO::NGSPipeline::RNAseq::Common
            CO::NGSPipeline::Common);

sub new {
	my $class = shift;
	$class = ref($class) ? ref($class) : $class;
	
	my $pipeline = shift;
	
	my $self = {"pipeline" => $pipeline};
	
	return bless $self, $class;
}

sub align {
	my $self = shift;
	
	my %param = ( "fastq1" => undef,
	              "fastq2" => undef,
				  "sample_id" => "sample",
				  "delete_input" => 0,
				  "strand" => 0,
				  @_);
	
	my $fastq1 = to_abs_path($param{fastq1});
	my $fastq2 = to_abs_path($param{fastq2});
	my $output = to_abs_path($param{output});
	my $delete_input = $param{delete_input};
	my $sample_id = $param{sample_id};
	my $strand = $param{strand};
	
	unless($outptu =~/\.bam$/) {
		die "Only permit outputting bam file in STAR alignment.\n";
	}
	
	my $pipeline = $self->{pipeline};
	
	$pipeline->add_command("mkfifo $pipeline->{dir}/fastq1 $pipeline->{dir}/fastq2", 0);
    $pipeline->add_command("zcat -c $fastq1 > $pipeline->{dir}/fastq1 &", 0);
    $pipeline->add_command("zcat -c $fastq2 > $pipeline->{dir}/fastq2 &", 0);
	
	if($strand) {
		$pipeline->add_command("STAR-2.3.0e --genomeDir $STAR_GENOME --readFilesIn $pipeline->{dir}/fastq1 $pipeline->{dir}/fastq2 --runThreadN 8 --outFileNamePrefix $pipeline->{dir}/$sample_id. --genomeLoad LoadAndRemove --alignIntronMax 500000 --alignMatesGapMax 500000 --outSAMunmapped Within --outFilterMultimapNmax 2 --outFilterMismatchNmax 3 --outFilterMismatchNoverLmax 0.3 --sjdbOverhang 50 --chimSegmentMin 15 --chimScoreMin 1 --chimScoreJunctionNonGTAG 0 --chimJunctionOverhangMin 15 --outStd SAM | mbuffer -q -m 2G -l /dev/null | samtools view -uSbh - | mbuffer -q -m 2G -l /dev/null | samtools sort - $pipeline->{dir}/tmp_$sample_id");
	} else {
		$pipeline->add_command("STAR-2.3.0e --outSAMstrandField intronMotif --genomeDir $STAR_GENOME --readFilesIn $pipeline->{dir}/fastq1 $pipeline->{dir}/fastq2 --runThreadN 8 --outFileNamePrefix $pipeline->{dir}/$sample_id. --genomeLoad LoadAndRemove --alignIntronMax 500000 --alignMatesGapMax 500000 --outSAMunmapped Within --outFilterMultimapNmax 2 --outFilterMismatchNmax 3 --outFilterMismatchNoverLmax 0.3 --sjdbOverhang 50 --chimSegmentMin 15 --chimScoreMin 1 --chimScoreJunctionNonGTAG 0 --chimJunctionOverhangMin 15 --outStd SAM | mbuffer -q -m 2G -l /dev/null | samtools view -uSbh - | mbuffer -q -m 2G -l /dev/null | samtools sort - $pipeline->{dir}/tmp_$sample_id");	
	}
	$pipeline->add_command("rm $pipeline->{dir}/fastq1", 0);
    $pipeline->add_command("rm $pipeline->{dir}/fastq2", 0);
	$pipeline->add_command("mv $pipeline->{dir}/tmp_$sample_id.bam $output", 0);
	
	$pipeline->check_filesize("$output");
	my $qid = $pipeline->run("-N" => $pipeline->get_job_name ? $pipeline->get_job_name : "_star_align",
							 "-l" => { nodes => "1:ppn=8:lsdf", 
									    mem => "40GB",
										walltime => "10:00:00"});

	return($qid);

}

1;
