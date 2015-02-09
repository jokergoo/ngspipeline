use strict;

my @files = @ARGV;

foreach my $file (@files) {
	open PIPE, "zcat $file |";
	while(<PIPE>) {
		print $_;
	}
	close PIPE;
}