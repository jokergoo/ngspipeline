#######################################
#
# Collect fastqc results
#
# Zuguang Gu
#
#######################################

use strict;
use File::Basename;
use Digest::MD5 qw(md5 md5_hex md5_base64);
use File::Spec;
use Getopt::Long;

if(scalar(@ARGV) == 0) {
print STDERR "examples:
perl collect_fastqc.pl --path \"/icgc/dkfzlsdf/analysis/hipo/hipo_047/whole_genome_sequencing/*/fastqc/*\" \\
                       --regexp \"\\/fastqc\\/(.*?_WHOLE_GENOME_PAIRED)\"
perl collect_fastqc.pl --path \"/icgc/dkfzlsdf/analysis/hipo/hipo_047/whole_genome_sequencing/*/fastqc/*\" \\ 
                       --regexp \"\\/fastq(.*?)_(WHOLE_GENOME)_PAIRED\" \\
                       --capture \"\\\$1-\\\$2\" 
perl collect_fastqc.pl --path \"/icgc/dkfzlsdf/analysis/hipo/hipo_016/chipseq_gbm/*/fastqc*/*\" \\
                       --regexp \"chipseq_gbm\\/(.*fastqc.*?)\\/\" \\
                       --output fastqc_report
";
exit 123;
}


my $output_dir;
my $path;
my $regexp;
my $capture;
GetOptions("output=s" => \$output_dir,
	       "path=s" => \$path,
	       "regexp=s" => \$regexp,
	       "capture=s" => \$capture) or die "opt wrong\n";


if(!defined($output_dir)) {
	$output_dir = "fastqc_collect";
}
if(!defined($capture)) {
	$capture = "\$1";
}

$output_dir = to_abs_path($output_dir);

`mkdir $output_dir` if(!-e $output_dir);

open HTML, ">$output_dir/index.html";

print HTML "<html>
<head>
<style>
table {
	border-top: 1px solid black;
	border-left: 1px solid black;
}

td, th {
	border-right: 1px solid black;
	border-bottom: 1px solid black;
}

.PASS, .WARN, .FAIL {
	text-align:center;
}

.PASS {
	color:green;
}

.WARN {
	color:orange;
}

.FAIL {
	color:red;
}
</style>
</head>
<body>
";

print HTML "<table cellspacing=0 cellpadding=0>\n";
print HTML "<tr><th>&nbsp;</th><th>Basic Statistics</th>
<th>Per base sequence quality</th>
<th>Per sequence quality scores</th>
<th>Per base sequence content</th>
<th>Per base GC content </th>
<th>Per sequence GC content </th>
<th>Per base N content</th>
<th>Sequence Length Distribution</th>
<th>Sequence Duplication Levels</th>
<th>Overrepresented sequences</th>
<th>Kmer Content</th>
</tr>\n";

my $icons = {
	PASS => "tick.png",
	WARN => "warning.png",
	FAIL => "error.png"
};

my $cmd =  "ls -d $path";
my @folder = `$cmd`;

@folder = map {$_=~s/\n//msg;$_} @folder;
@folder = grep {-d} @folder;


if(! scalar(@folder)) {
	die "no directory was found.\n";
}

my $pid_eval = "\"$capture\"";


for my $d (@folder) {
	my $regexp = qr/$regexp/;
	if($d =~/$regexp/) {
		my $sample_name = eval($pid_eval);
		my $new_d = md5_hex($d);
		print STDERR "copy $d to $new_d\n";
		`cp -r $d $output_dir/$new_d`;
		print_tr($sample_name, "$output_dir/$new_d", $new_d);
	} else {
		die "cannot find sample name ($d)\n";
	}
}


print HTML "</table>\n";
print HTML "</body>\n</html>\n";
close HTML;
`tar -zcvf $output_dir.tar.gz $output_dir`;

sub print_tr {
	my $sample_name = shift;
	my $fastqc_dir = shift;
	my $hash = shift;
	
	print HTML "<tr><td><a href='$hash/fastqc_report.html'>$sample_name</a></td>";
	open F, "$fastqc_dir/summary.txt";
	my $i_line = 0;
	while(my $line = <F>) {
		my @tmp = split "\t", $line;
		print HTML "<td class='$tmp[0]'><a href='$hash/fastqc_report.html#M$i_line'><img src='$hash/Icons/$icons->{$tmp[0]}' alt='$tmp[0]' /></a></td>";
		$i_line ++;
	}
	print HTML "</tr>\n";
}

sub to_abs_path {
	my $path = shift;
	if(defined($path)) {
		return File::Spec->file_name_is_absolute($path) ? $path : File::Spec->rel2abs($path);
	} else {
		return undef;
	}
}
