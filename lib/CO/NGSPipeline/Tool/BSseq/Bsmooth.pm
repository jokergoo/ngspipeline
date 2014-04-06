package CO::NGSPipeline::Tool::BSseq::Bsmooth;

use strict;
use CO::NGSPipeline::Tool::BSseq::Config;
use CO::NGSPipeline::Tool::Config;
use CO::Utils;
use File::Temp qw/tempfile/;
use File::Basename;

use base qw(CO::NGSPipeline::Tool::BSseq::Common
            CO::NGSPipeline::Tool::Common
			CO::NGSPipeline);

sub new {
	my $class = shift;
	$class = ref($class) ? ref($class) : $class;
	
	my $self = {};
	
	return bless $self, $class;
}


sub RData {
	my $self = shift;
	
	my %param = ( "input" => undef,
		          "chr" => undef,
	              "sample_id" => undef,
				  "output" => undef,
				  "delete_input" => 0,
				  @_);

	my $input = $param{input};
	my $sampleName = $param{sample_id};
	my $delete_input = $param{delete_input};
	my $chr = $param{chr};
	my $output = to_abs_path($param{output});
	
	my $output2 = $output;
	$output2 =~s/RData$/smoothed.RData/i;
	
	my $pm = $self->get_pipeline_maker;
	
	my $input_chr = "$pm->{dir}/".basename($input).".$chr";
	
	$pm->add_command("grep -e \"^$chr\\s\" $input > $input_chr", 0);
	
	my $file_str = "c(\"$input_chr\")";
	my $sampleName_str = "c(\"$sampleName\")";
	
	my $R = <<R;
	
library(bsseq, lib.loc = "/home/guz/R/x86_64-unknown-linux-gnu-library/3.0")
d = read.table($file_str, sep = "\t", stringsAsFactors = FALSE)
chr = d[[1]]
M = cbind(round(d[[4]]/100*d[[5]]))
Cov = cbind(d[[5]])
bsseq_$chr = BSseq(chr = chr, pos = d[[2]], M = M, Cov = Cov, sampleNames = $sampleName_str)

save(bsseq_$chr, file = "$output")

smooth_$chr = BSmooth(bsseq_$chr)

save(smooth_$chr, file = "$output2")

R

	
	
	my ($tmp_fh, $tmp_file) = tempfile(DIR => $pm->{tmp_dir}, SUFFIX => '.R');
	
	$tmp_file = to_abs_path($tmp_file);
	
	print $tmp_fh $R;
	close $tmp_fh;
	
	$pm->add_command("Rscript $tmp_file");
	$pm->del_file($input_chr);
	$pm->del_file($tmp_file);
	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_dmr_extract",
							 "-l" => { nodes => "1:ppn=1:lsdf", 
									    mem => "30GB",
										walltime => "30:00:00"});

}

1;
