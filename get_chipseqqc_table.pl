use strict;
use File::Basename;

my @f = `ls /icgc/dkfzlsdf/analysis/hipo/hipo_016/chipseq_gbm/*/*mkdup.flagstat`;

print "PID\ttotal_reads\tmapping_rate\tlanes\tduplication_rate\n";
for my $f (@f) {

	my $dir = dirname($f);
	#print "ls $dir/*.flagstat";
	my @tmp = `ls $dir/*.flagstat`;
	my $n_lanes = scalar(@tmp) - 1;
	open F, $f;
	
	$f =~/chipseq_gbm\/(.*?)\//;
	my $pid = $1;
	
	print $pid, "\t";
	
	my $line = <F>;
	$line =~/^(\d+)/;
	print $1, "\t";
	
	$line = <F>;
	$line = <F>;
	$line =~/(\d+\.\d+%)/;
	print $1, "\t";
	print $n_lanes, "\t";
	open F, "$dir/$pid.mkdup.mkdup.metrics" or die "cannot open $dir/$pid.mkdup.mkdup.metrics";
	while($line = <F>) {
	 last if $line =~/^LIBRARY/;
	}
	
	$line = <F>;
	@tmp = split " ", $line;
	print $tmp[8], "\n";
}