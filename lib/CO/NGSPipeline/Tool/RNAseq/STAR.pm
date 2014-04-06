package CO::NGSPipeline::Tool::RNAseq::STAR;

use strict;
use CO::NGSPipeline::Tool::RNAseq::Config;
use CO::NGSPipeline::Tool::Config;
use CO::Utils;

use base qw(CO::NGSPipeline::Tool::RNAseq::Common
            CO::NGSPipeline::Tool::Common
			CO::NGSPipeline);

sub new {
	my $class = shift;
	$class = ref($class) ? ref($class) : $class;
	
	my $self = {};
	
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
	
	unless($output =~/\.bam$/) {
		die "Only permit outputting bam file in STAR alignment.\n";
	}
	
	my $pm = $self->get_pipeline_maker;
	
	my $rt = int(rand(999999)).time();
	$pm->add_command("mkfifo $pm->{dir}/$rt.fastq1 $pm->{dir}/$rt.fastq2", 0);
	 $pm->add_command("zcat -c $fastq1 > $pm->{dir}/$rt.fastq1 &", 0);
    $pm->add_command("zcat -c $fastq2 > $pm->{dir}/$rt.fastq2 &", 0);
	
	my $r = time().rand();
	if($strand) {
		$pm->add_command("$STAR --genomeDir $STAR_GENOME --readFilesIn $pm->{dir}/$rt.fastq1 $pm->{dir}/$rt.fastq2 --runThreadN 8 --outFileNamePrefix $pm->{dir}/$sample_id. --genomeLoad LoadAndRemove --alignIntronMax 500000 --alignMatesGapMax 500000 --outSAMunmapped Within --outFilterMultimapNmax 2 --outFilterMismatchNmax 3 --outFilterMismatchNoverLmax 0.3 --sjdbOverhang 50 --chimSegmentMin 15 --chimScoreMin 1 --chimScoreJunctionNonGTAG 0 --chimJunctionOverhangMin 15 --outStd SAM | mbuffer -q -m 2G -l /dev/null | $SAMTOOLS view -uSbh - | mbuffer -q -m 2G -l /dev/null | $SAMTOOLS sort - $pm->{dir}/$r");
	} else {
		$pm->add_command("$STAR --outSAMstrandField intronMotif --genomeDir $STAR_GENOME --readFilesIn $pm->{dir}/$rt.fastq1 $pm->{dir}/$rt.fastq2 --runThreadN 8 --outFileNamePrefix $pm->{dir}/$sample_id. --genomeLoad LoadAndRemove --alignIntronMax 500000 --alignMatesGapMax 500000 --outSAMunmapped Within --outFilterMultimapNmax 2 --outFilterMismatchNmax 3 --outFilterMismatchNoverLmax 0.3 --sjdbOverhang 50 --chimSegmentMin 15 --chimScoreMin 1 --chimScoreJunctionNonGTAG 0 --chimJunctionOverhangMin 15 --outStd SAM | mbuffer -q -m 2G -l /dev/null | $SAMTOOLS view -uSbh - | mbuffer -q -m 2G -l /dev/null | $SAMTOOLS sort - $pm->{dir}/$r");	
	}
	$pm->add_command("rm $pm->{dir}/$rt.fastq1", 0);
    $pm->add_command("rm $pm->{dir}/$rt.fastq2", 0);
	$pm->add_command("mv $pm->{dir}/$r.bam $output", 0);
	$pm->add_command("$SAMTOOLS index $output");
	
	$pm->check_filesize("$output");
	$pm->del_file($fastq1, $fastq2) if($delete_input);
	
	
	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_star_align",
							 "-l" => { nodes => "1:ppn=8:lsdf", 
									    mem => "40GB",
										walltime => "10:00:00"});

	return($qid);

}

1;
