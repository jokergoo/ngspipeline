package CO::NGSPipeline::BSseq::RData;

use strict;
use CO::NGSPipeline::BSseq::Config;
use CO::NGSPipeline::Utils;
use File::Temp qw/tempfile/;
use File::Basename;

use base qw(CO::NGSPipeline::BSseq::Common
            CO::NGSPipeline::Common);

sub new {
	my $class = shift;
	$class = ref($class) ? ref($class) : $class;
	
	my $pipeline = shift;
	
	my $self = {"pipeline" => $pipeline};
	
	return bless $self, $class;
}


sub RData {
	my $self = shift;
	
	my %param = ( "input" => undef,
		          "chr" => undef,
	              "sampleName" => undef,
				  "sampleClass" => undef,
				  "output" => undef,
				  "delete_input" => 0,
				  @_);

	my $input = $param{input};
	my $sampleName = $param{sampleName};
	my $sampleClass = $param{sampleClass};
	my $delete_input = $param{delete_input};
	my $chr = $param{chr};
	my $output = to_abs_path($param{output});
	
	my $pipeline = $self->{pipeline};
	
	my $input_chr = "$pipeline->{dir}/".basename($input).".$chr";
	
	$pipeline->add_command("grep -e \"^$chr\\s\" $input > $input_chr", 0);
	
	my $file_str = "c(\"$input_chr\")";
	my $sampleName_str = "c(\"$sampleName\")";
	my $sampleClass_str = "c(\"$sampleClass\")";
	
	my $R = <<R;
	
library(bsseq)
#bsseq_$chr = read.bismark(files = $file_str, sampleNames = $sampleName_str, rmZeroCov = TRUE, verbose=TRUE)
d = read.table($file_str, sep = "\t", stringsAsFactors = FALSE)
chr = d[[1]]
M = cbind(round(d[[4]]/100*d[[5]]))
Cov = cbind(d[[5]])
bsseq_$chr = BSseq(chr = chr, pos = d[[2]], M = M, Cov = Cov, sampleNames = $sampleName_str)

save(bsseq_$chr, file = "$output")

R

	
	
	my ($tmp_fh, $tmp_file) = tempfile(DIR => $pipeline->{tmp_dir}, SUFFIX => '.R');
	
	$tmp_file = to_abs_path($tmp_file);
	
	print $tmp_fh $R;
	close $tmp_fh;
	
	$pipeline->add_command("Rscript $tmp_file");
	$pipeline->del_file($input_chr);
	$pipeline->del_file($tmp_file);
	my $qid = $pipeline->run("-N" => $pipeline->get_job_name ? $pipeline->get_job_name : "_dmr_extract",
							 "-l" => { nodes => "1:ppn=1:lsdf", 
									    mem => "30GB",
										walltime => "10:00:00"});

}

1;
