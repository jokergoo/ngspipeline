package CO::NGSPipeline::Tool::BSseq::BSMAP;

use strict;
use CO::NGSPipeline::Tool::BSseq::Config;
use CO::NGSPipeline::Tool::Config;
use CO::Utils;

use base qw(CO::NGSPipeline::Tool::BSseq::Common
            CO::NGSPipeline::Tool::Common
			CO::NGSPipeline);

sub new {
	my $class = shift;
	$class = ref($class) ? ref($class) : $class;
	
	my $self = {};
	
	return bless $self, $class;
}

=item C<$self-E<gt>align(HASH)>

BSMAP alignment

  fastq1        path of read 1
  fastq2        path of read 2
  output        path of output
  delete_input  whether delete input files
 
=cut
sub align {
	my $self = shift;
	
	my %param = ( "fastq1" => undef,
	              "fastq2" => undef,
				  "output" => undef,
				  "delete_input" => 0,
				  "species" => 'human',
				  @_);
	
	my $fastq1 = to_abs_path($param{fastq1});
	my $fastq2 = to_abs_path($param{fastq2});
	my $output = to_abs_path($param{output});
	my $delete_input = $param{delete_input};
	my $species = $param{species};

	if(lc($species) eq "mouse") {
		$BSMAP_REF_GENOME = $BSMAP_REF_GENOME_MM;
	}
	
	unless($output =~/\.bam$/) {
		die "Only permit outputting bam file in BSMAP alignment.\n";
	}
	
	my $pm = $self->get_pipeline_maker;
	
	$pm->add_command("PATH=\$PATH:$BSMAP_BIN_DIR/", 0);
	my $random_str = time().int(rand(9999999));
	$pm->add_command("mkfifo \$PBS_SCRATCH_DIR/\$PBS_JOBID/$random_str.sam", 0);

	# note: when this pipe fails, it will not return a failed flag, only exit!
	$pm->add_command("bsmap -a $fastq1 -b $fastq2 -d $BSMAP_GENOME_DIR/$BSMAP_REF_GENOME -o \$PBS_SCRATCH_DIR/\$PBS_JOBID/$random_str.sam -p 8 -r 0 -v 8 &"); # unique match, at most 8 mismatches
	
	my $output_base = $output; $output_base =~s/\.bam$//;
	$pm->add_command("$SAMTOOLS view -uSbh \$PBS_SCRATCH_DIR/\$PBS_JOBID/$random_str.sam | mbuffer -q -m 2G -l /dev/null | $SAMTOOLS sort - $output_base");
	$pm->add_command("$SAMTOOLS index $output");
	$pm->add_command("rm \$PBS_SCRATCH_DIR/\$PBS_JOBID/$random_str.sam", 0);
	
	
	$pm->check_filesize($output); # 1M
	$pm->del_file($fastq1, $fastq2) if($delete_input);
	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_bsmap_align",
							 "-l" => { nodes => "1:ppn=8:lsdf", 
									    mem => "12GB",
										walltime => "150:00:00"});

	return($qid);
}

sub lambda_conversion {
	my $self = shift;
	
	my %param = ( "bam" => undef,
	               "output" => undef,
	               "delete_input" => 0,
	               "species" => 'human',
				   @_);
	
	my $bam    = to_abs_path($param{bam});
	my $output = to_abs_path($param{output});
	my $delete_input = $param{delete_input};
	my $species = $param{species};

	if(lc($species) eq "mouse") {
		$BSMAP_REF_GENOME = $BSMAP_REF_GENOME_MM;
	}
	
	my $pm = $self->get_pipeline_maker;
	
	$pm->add_command("$SAMTOOLS view -h $bam lambda > $pm->{tmp_dir}/lambda.sam");
	$pm->add_command("python $BSMAP_BIN_DIR/methratio.py -d $BSMAP_GENOME_DIR/$BSMAP_REF_GENOME --chr=lambda -o $pm->{tmp_dir}/lambda.out -p $pm->{tmp_dir}/lambda.sam -z");
	$pm->add_command("cat $pm->{tmp_dir}/lambda.out | awk '{a+=\$5}END{print 1-a/(NR-1)}' > $output", 0);
	
	$pm->del_file("$pm->{tmp_dir}/lambda.sam");
	$pm->del_file("$pm->{tmp_dir}/lambda.out");

	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_bsmap_lambda_conversion",
							 "-l" => { nodes => "1:ppn=1:lsdf", 
									    mem => "1GB",
										walltime => "20:00:00"});

	return($qid);

}

=item C<$self-E<gt>call_methylation(HASH)>

BSMAP methylation calling

  sam           sam files
  output        output
  delete_input  whether delete input files
 
=cut
sub call_methylation {
	my $self = shift;
	
	my %param = ( "sam" => undef,
	               "output" => undef,
	               "delete_input" => 0,
	               "species" => 'human',
				   @_);
	
	my $sam    = to_abs_path($param{sam});
	my $output = to_abs_path($param{output});
	my $delete_input = $param{delete_input};
	my $species = $param{species};

	if(lc($species) eq "mouse") {
		$BSMAP_REF_GENOME = $BSMAP_REF_GENOME_MM;
	}
	
	my $pm = $self->get_pipeline_maker;
	
	$pm->add_command("python $BSMAP_BIN_DIR/methratio.py -d $BSMAP_GENOME_DIR/$BSMAP_REF_GENOME -o $output -p $sam");

	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_bsmap_call_methylation",
							 "-l" => { nodes => "1:ppn=1:lsdf", 
									    mem => "40GB",
										walltime => "100:00:00"});

	return($qid);
}

1;
