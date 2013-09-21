package CO::NGSPipeline::Tool::BSseq::methylCtools;

use strict;
use CO::NGSPipeline::Tool::BSseq::Config;
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




sub fqconv {
	my $self = shift;
	
	my %param = ( "fastq" => undef,
	              "which_pair" => 1,
	              "output" => undef,
				  "delete_input" => 0,
				  @_);
	
	my $fastq = to_abs_path($param{fastq});
	my $output = to_abs_path($param{output});
	my $which_pair = $param{which_pair};
	my $delete_input = $param{delete_input};
	
	my $pm = $self->get_pipeline_maker;
	
	if($which_pair == 1) {
		$pm->add_command("$METHYLCTOOLS_BIN_DIR/methylCtools fqconv -1 $fastq $output");
	} else {
		$pm->add_command("$METHYLCTOOLS_BIN_DIR/methylCtools fqconv -2 $fastq $output");
	}
	$pm->del_file($fastq) if($delete_input);
	$pm->check_filesize($output); # 1M
	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_methylCtools_fqconv",
							 "-l" => { nodes => "1:ppn=1:lsdf", 
									    mem => "10GB",
										walltime => "20:00:00"});

	return($qid);
}



sub bconv {
	my $self = shift;
	
	my %param = ( "bam" => undef,
	              "output" => undef,
				  "delete_input" => 0,
				  @_);
	
	my $bam = to_abs_path($param{bam});
	my $output = to_abs_path($param{output});
	my $delete_input = $param{delete_input};
	
	my $pm = $self->get_pipeline_maker;
	
	$pm->add_command("$METHYLCTOOLS_BIN_DIR/methylCtools bconv $bam $output");
	$pm->del_file($bam) if($delete_input);
	$pm->check_filesize($output); # 1M
	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_methylCtools_bconv",
							 "-l" => { nodes => "1:ppn=1:lsdf", 
									    mem => "10GB",
										walltime => "50:00:00"});

	return($qid);
}


sub bcall {
	my $self = shift;
	
	my %param = ( "bam" => undef,
	              "output" => undef,
				  "delete_input" => 0,
				  @_);
	
	my $bam = to_abs_path($param{bam});
	my $output = to_abs_path($param{output});
	my $delete_input = $param{delete_input};
	
	my $pm = $self->get_pipeline_maker;
	
	my $bai = $bam;
	$bai =~s/\.bam$/.bai/;
	
	$pm->add_command("ln -s $bai $bam.bai", 0);
 	$pm->add_command("$METHYLCTOOLS_BIN_DIR/methylCtools bcall -t -m $output.metrics.txt $METHYLCTOOLS_GENOME_DIR/$METHYLCTOOLS_REFERENCE_POS_FILE $bam - | grep -v 'lambda' | bgzip > $output");
	$pm->add_command("tabix -s 1 b 2 -e 2 $output");
	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_methylCtools_bcall",
							 "-l" => { nodes => "1:ppn=1:lsdf", 
									    mem => "10GB",
										walltime => "150:00:00"});

	return($qid);
}

sub lambda_conversion {
	my $self = shift;
	
	my %param = ( "bam" => undef,
	               "output" => undef,
	               "delete_input" => 0,
				   @_);
	
	my $bam    = to_abs_path($param{bam});
	my $output = to_abs_path($param{output});
	my $delete_input = $param{delete_input};
	
	my $pm = $self->get_pipeline_maker;
	
	$pm->add_command("samtools view -hb $bam lambda > $pm->{tmp_dir}/lambda.bam");
	$pm->add_command("samtools index $pm->{tmp_dir}/lambda.bam");
	$pm->add_command("$METHYLCTOOLS_BIN_DIR/methylCtools bcall -r lambda $METHYLCTOOLS_GENOME_DIR/$METHYLCTOOLS_REFERENCE_POS_FILE $pm->{tmp_dir}/lambda.bam - | bgzip > $pm->{tmp_dir}/lambda.bcall.gz");
	$pm->add_command("zcat $pm->{tmp_dir}/lambda.bcall.gz | awk '{a+=\$5/(\$5+\$6)}END{print 1-a/NR}' > $output", 0);
	
	$pm->del_file("$pm->{tmp_dir}/lambda.bam");
	$pm->del_file("$pm->{tmp_dir}/lambda.bam.bai");
	$pm->del_file("$pm->{tmp_dir}/lambda.out");

	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_bsmap_lambda_conversion",
							 "-l" => { nodes => "1:ppn=1:lsdf", 
									    mem => "1GB",
										walltime => "20:00:00"});

	return($qid);

}


1;
