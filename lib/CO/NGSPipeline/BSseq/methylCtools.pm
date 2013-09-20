package CO::NGSPipeline::BSseq::methylCtools;

use strict;
use CO::NGSPipeline::BSseq::Config;
use CO::NGSPipeline::Utils;

use base qw(CO::NGSPipeline::BSseq::Common
            CO::NGSPipeline::Common);

sub new {
	my $class = shift;
	$class = ref($class) ? ref($class) : $class;
	
	my $pipeline = shift;
	
	my $self = {"pipeline" => $pipeline};
	
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
	
	my $pipeline = $self->{pipeline};
	
	if($which_pair == 1) {
		$pipeline->add_command("$METHYLCTOOLS_BIN_DIR/methylCtools fqconv -1 $fastq $output");
	} else {
		$pipeline->add_command("$METHYLCTOOLS_BIN_DIR/methylCtools fqconv -2 $fastq $output");
	}
	$pipeline->del_file($fastq) if($delete_input);
	$pipeline->check_filesize($output); # 1M
	my $qid = $pipeline->run("-N" => $pipeline->get_job_name ? $pipeline->get_job_name : "_methylCtools_fqconv",
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
	
	my $pipeline = $self->{pipeline};
	
	$pipeline->add_command("$METHYLCTOOLS_BIN_DIR/methylCtools bconv $bam $output");
	$pipeline->del_file($bam) if($delete_input);
	$pipeline->check_filesize($output); # 1M
	my $qid = $pipeline->run("-N" => $pipeline->get_job_name ? $pipeline->get_job_name : "_methylCtools_bconv",
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
	
	my $pipeline = $self->{pipeline};
	
	my $bai = $bam;
	$bai =~s/\.bam$/.bai/;
	
	$pipeline->add_command("ln -s $bai $bam.bai", 0);
 	$pipeline->add_command("$METHYLCTOOLS_BIN_DIR/methylCtools bcall -t -m $output.metrics.txt $METHYLCTOOLS_GENOME_DIR/$METHYLCTOOLS_REFERENCE_POS_FILE $bam - | grep -v 'lambda' | bgzip > $output");
	$pipeline->add_command("tabix -s 1 b 2 -e 2 $output");
	my $qid = $pipeline->run("-N" => $pipeline->get_job_name ? $pipeline->get_job_name : "_methylCtools_bcall",
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
	
	my $pipeline = $self->{pipeline};
	
	$pipeline->add_command("samtools view -hb $bam lambda > $pipeline->{tmp_dir}/lambda.bam");
	$pipeline->add_command("samtools index $pipeline->{tmp_dir}/lambda.bam");
	$pipeline->add_command("$METHYLCTOOLS_BIN_DIR/methylCtools bcall -r lambda $METHYLCTOOLS_GENOME_DIR/$METHYLCTOOLS_REFERENCE_POS_FILE $pipeline->{tmp_dir}/lambda.bam - | bgzip > $pipeline->{tmp_dir}/lambda.bcall.gz");
	$pipeline->add_command("zcat $pipeline->{tmp_dir}/lambda.bcall.gz | awk '{a+=\$5/(\$5+\$6)}END{print 1-a/NR}' > $output", 0);
	
	$pipeline->del_file("$pipeline->{tmp_dir}/lambda.bam");
	$pipeline->del_file("$pipeline->{tmp_dir}/lambda.bam.bai");
	$pipeline->del_file("$pipeline->{tmp_dir}/lambda.out");

	my $qid = $pipeline->run("-N" => $pipeline->get_job_name ? $pipeline->get_job_name : "_bsmap_lambda_conversion",
							 "-l" => { nodes => "1:ppn=1:lsdf", 
									    mem => "1GB",
										walltime => "20:00:00"});

	return($qid);

}


1;
