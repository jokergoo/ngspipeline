package CO::NGSPipeline::Tool::RNAseq::TopHat;

use strict;
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
				  "output" => undef,
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
	
	my $pm = $self->get_pipeline_maker;
	
	if($strand) {
		$pm->add_command("tophat2 -o $pm->{dir} -p 8 --library-type fr-firststrand -r 200 --mate-std-dev 50 --b2-sensitive -g 1 --no-coverage-search --GTF $GENCODE_GTF --transcriptome-index=$GENCODE_BOWTIE2_INDEX $BOWTIE2_INDEX $fastq1 $fastq2");
	} else {
		$pm->add_command("tophat2 -o $pm->{dir} -p 8 --library-type fr-unstranded -r 200 --mate-std-dev 50 --b2-sensitive -g 1 --no-coverage-search --GTF $GENCODE_GTF --transcriptome-index=$GENCODE_BOWTIE2_INDEX $BOWTIE2_INDEX $fastq1 $fastq2");
	}
	$pm->add_command("mv $pm->{dir}/accepted_hits.bam $output", 0);
	$pm->check_filesize("$output");
	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_tophat2_align",
							 "-l" => { nodes => "1:ppn=8:lsdf", 
									    mem => "10GB",
										walltime => "150:00:00"});

	return($qid);

}

1;
