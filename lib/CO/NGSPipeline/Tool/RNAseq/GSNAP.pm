package CO::NGSPipeline::Tool::RNAseq::GSNAP;

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
	
	unless($output =~/\.bam$/) {
		die "Only permit outputting bam file in GSNAP alignment.\n";
	}
	
	my $pm = $self->get_pipeline_maker;
	
	# which value should be set to -B to ensure maximum memory usage under 90G
	$pm->add_command("gsnap -D $GSNAP_GENOME_DIR -d $GSNAP_GENOME --nthreads=16 -B 5 -s $GSNAP_IIT -n 2 -Q --nofails --format=sam --gunzip $fastq1 $fastq2 | mbuffer -q -m 2G -l /dev/null | samtools view -uSbh - | mbuffer -q -m 2G -l /dev/null | samtools sort - $pm->{dir}/tmp_$sample_id");

	$pm->add_command("mv $pm->{dir}/tmp_$sample_id.bam $output", 0);
	$pm->check_filesize("$output");
	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_gsnap_align",
							 "-l" => { nodes => "1:ppn=16:lsdf", 
									    mem => "20GB",
										walltime => "300:00:00"});

	return($qid);

}

1;
