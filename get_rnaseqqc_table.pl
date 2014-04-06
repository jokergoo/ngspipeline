use strict;

my $dir = $ARGV[0];

chdir $dir;

my @sample = glob("*");

my $h = {};
print "\tTotal\tUnique\tDuplicates\tDuplication Rate\tMapped\tMapping Rate\tMapped Unique\tMapped Unique Rate\trRNA rate\tIntragenic Rate\tExonic Rate\tIntronic Rate\tIntergenic Rate\tEnd 1 % Sense\tEnd 2 % Sense\n";
foreach my $sid (sort @sample) {

	next if($sid =~/old/);

	my $d = read_data($sid);

	print "$sid\t";
	print join "\t", @$d;
	print "\n";
}




sub read_data {
	my $sid = shift;

	my $path = "$dir/$sid/rnaseqqc/$sid/$sid.metrics.txt";
	open F, $path or die "cannot open $path\n";

	my @d;

	<F>;
	my $line = <F>; chomp $line;
	my @tmp = split "\t", $line;

	push(@d, @tmp[0..3]);

	<F>;
	$line = <F>; chomp $line;
	@tmp = split "\t", $line;

	push(@d, @tmp[0..3]); push(@d, $tmp[5]);

	<F>;
	$line = <F>; chomp $line;
	@tmp = split "\t", $line;

	push(@d, @tmp[0..3]);

	<F>;
	$line = <F>; chomp $line;
	@tmp = split "\t", $line;

	push(@d, @tmp[4..5]);

	return \@d;
}