package CO::NGSPipeline::Tool::ChIPSeq::Common;

##############################################################################
# provide command method for each pipeline. It is a base module for all specific
# pipelines.

use strict;
use CO::NGSPipeline::Tool::ChIPSeq::Config;
use CO::NGSPipeline::Tool::Config;
use CO::Utils;
use File::Basename;

use base qw(CO::NGSPipeline::Tool::Common
			CO::NGSPipeline);

sub new {
	my $class = shift;
	$class = ref($class) ? ref($class) : $class;

	my $self = {};
	
	return bless $self, $class;
}

sub data_transform {
	my $self = shift;
	
	my %param = ( "sam" => undef,
				  @_);
	
	my $sam    = to_abs_path( $param{sam} );
	my $bed = $sam; $bed =~s/\.(sam|bam)$/.bed/;
	my $tdf = $sam; $tdf =~s/\.(sam|bam)$/.tdf/;
	
	my $pm = $self->get_pipeline_maker;
	
	$pm->add_command("bedtools bamtobed -i $sam > $bed");
	$pm->add_command("igvtools count --minMapQuality 1 $sam $tdf hg19");
	
	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_data_transform",
							 "-l" => { nodes => "1:ppn=1:lsdf", 
									    mem => "1GB",
										walltime => "1:00:00"});
	return($qid);
}

sub peak_calling {
	my $self = shift;

	my %param = ("sam" => undef,
				"output_dir" => undef,
				"sample_id" => undef,
				"homer_qc_dir" => undef,
		          @_);

	my $sam    = to_abs_path( $param{sam} );
	my $sample_id = $param{sample_id};
	my $output_dir = to_abs_path( $param{output_dir} );
	my $homer_qc_dir = to_abs_path( $param{homer_qc_dir} );
	my $pm = $self->get_pipeline_maker;
	
	$pm->add_command("fragment_length=`cat $homer_qc_dir/tagInfo.txt | grep 'fragmentLengthEstimate' | cut -f2 -d=`", 0);
	$pm->add_command("half_fragment_length=`echo \$fragment_length/2 | bc`", 0);
	$pm->add_command("echo fragment_length=\$fragment_length", 0);
	$pm->add_command("echo half_fragment_length=\$half_fragment_length", 0);
	
	$pm->add_command("mkdir $output_dir/macs", 0);
	$pm->add_command("mkdir $output_dir/sicer", 0);
	$pm->add_command("cd $output_dir/macs", 0);
	$pm->add_command("rm -rf $sample_id"."_MACS_wiggle", 0);
	$pm->add_command("macs14 -t $sam -c /icgc/dkfzlsdf/analysis/hipo/hipo_016/chipseq_gbm/Pool_of_21_AK_samples/Pool_of_21_AK_samples.mkdup.bam --nomodel --nolambda --bw=\$fragment_length --tsize=51 --shiftsize=\$half_fragment_length --format=BAM --pvalue=1e-5 --name=$sample_id -g hs --wig");
	$pm->add_command("cd $output_dir/sicer; ln -s ../$sample_id.mkdup.bed; ln -s /icgc/dkfzlsdf/analysis/hipo/hipo_016/chipseq_gbm/Pool_of_21_AK_samples/Pool_of_21_AK_samples.mkdup.bed", 0);
	$pm->add_command("/ibios/co02/ishaque/prog/SICER_V1.1/SICER/SICER_nav_edit.v1.3.X.sh . $sample_id.mkdup.bed Pool_of_21_AK_samples.mkdup.bed . hg19 1 200 \$fragment_length 0.85");
	$pm->add_command("mkdir $output_dir/merged_sicer_macs", 0);
	$pm->add_command("cat $output_dir/macs/$sample_id"."_peaks.bed $output_dir/sicer/$sample_id.mkdup*W200-G200*0.05*broadPeak | cut -f 1,2,3 | sort -V -k1,1 -k2,2 | bedtools merge -i stdin > $output_dir/merged_sicer_macs/$sample_id.merged_sicer_macs.bed");

	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_peak_calling",
							 "-l" => { nodes => "1:ppn=1:lsdf", 
									    mem => "5GB",
										walltime => "50:00:00"});
	return($qid);
}



sub homer_qc {
	my $self = shift;

	my %param = ("sam" => undef,
				"output_dir" => undef,
				"sample_id" => undef,
		          @_);

	my $sam    = to_abs_path( $param{sam} );
	my $sample_id = $param{sample_id};
	my $output_dir = to_abs_path( $param{output_dir} );
	my $pm = $self->get_pipeline_maker;

	$pm->add_command("mkdir $output_dir/homer_qc", 0);
	$pm->add_command("/home/guz/soft/homer/bin/makeTagDirectory $output_dir/homer_qc $sam");

	my $cmd = "Rscript -e \"
pdf('$output_dir/homer_qc/$sample_id.homer_qc.pdf', width = 12, height = 12);
layout(rbind(1:2, c(3, 3)));
df = read.table('$output_dir/homer_qc/tagCountDistribution.txt', header = TRUE, sep = '\\t');
df = df[df[[1]]<=10 & df[[1]] > 0, ];
barplot(df[[2]], names.arg = df[[1]], col = 'blue', ylab = 'Fraction of total reads', xlab = 'Reads per position');


df = read.table('$output_dir/homer_qc/tagAutocorrelation.txt', header = TRUE, sep = '\\t');
plot(NULL, xlim = range(df[[1]]), ylim = range(c(df[[2]], df[[3]])), ylab = 'Total Read Paris', xlab = 'Relative Distance between Reads (bp)');
lines(df[[1]], df[[2]], col = 'blue');
lines(df[[1]], df[[3]], col = 'red');
legend('topright', lty = 1, col = c('blue', 'red'), legend = c('Same Strand', 'Opposite Strand'));
abline(h = 1, col = 'grey');
text(df[[1]][which.max(df[[3]])], max(df[[3]]), df[[1]][which.max(df[[3]])]);

source('~/project/development/meth_lib/common.R');
df = bedtools('bedtools genomecov -ibam $output_dir/$sample_id.mkdup.bam -g /ibios/co02/guz/MethSuite/bed_v19/hg19.len -bg');

df = df[df[[1]] %in% paste0('chr', c(1:22, 'X','Y')), ];
y = table(df[[4]]);
x = as.numeric(names(y));
plot(x, log10(as.vector(y)), type = 'h', ann = FALSE, axes = FALSE);
axis(side = 1);
axis(side = 2, at = 0:7, labels = sprintf('%.2e', 10^(0:7)/3095677412));
box();
title(xlab = 'per-base coverage', ylab = 'percentage in genome');


dev.off();
\"";
	$cmd =~s/\n//msg;
	$pm->add_command($cmd, 0);

	my $qid = $pm->run("-N" => $pm->get_job_name ? $pm->get_job_name : "_homer_qc",
							 "-l" => { nodes => "1:ppn=1:lsdf", 
									    mem => "2GB",
										walltime => "2:00:00"});
	return($qid);
}
1;
