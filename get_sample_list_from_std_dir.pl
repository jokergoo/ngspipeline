#!/bin/perl/bin/

BEGIN {
	use File::Basename;
	unshift(@INC, dirname($0)."/lib");
}

use strict;
use CO::NGSPipeline::Utils;
use File::Basename;
use List::Vectorize;

use Data::Dumper;

if(! @ARGV) {
	help_msg();
	exit 0;
}

my $list = read_std_dir($ARGV[0]);
my $pid = $list->{pid};
my $r1 = $list->{r1};
my $r2 = $list->{r2};

for(my $i = 0; $i < scalar(@{$pid}); $i ++) {
	for(my $j = 0; $j < scalar(@{$r1->[$i]}); $j ++) {
		print "$r1->[$i]->[$j]\t$r2->[$i]->[$j]\t$pid->[$i]\n";
	}

}

sub read_std_dir {
	my $dir = shift;
	
	if(! -d $dir) {
		help_msg();
		die "cannot find dir: $dir";
	}
	
	my $tree = dir_tree_as_hash($dir);
	
	if(!exists($tree->{"view-by-pid"})) {
		help_msg();
		die "cannot find view-by-pid directory in $dir\n";
	}
	
	my $pid = [ keys %{$tree->{"view-by-pid"}} ];
	my $r1 = [];
	my $r2 = [];
	
	for(my $i = 0; $i < scalar(@$pid); $i ++) {
		my $type = (keys %{$tree->{"view-by-pid"}->{$pid->[$i]}})[0];
		
		my @lanes = keys %{$tree->{"view-by-pid"}->{$pid->[$i]}->{$type}->{paired}};
		
		for(my $j = 0; $j < scalar(@lanes); $j ++) {
			($r1->[$i]->[$j], $r2->[$i]->[$j]) = keys %{$tree->{"view-by-pid"}->{$pid->[$i]}->{$type}->{paired}->{$lanes[$j]}->{sequence}};
		}
	}
	return {pid => $pid, r1 => $r1, r2 => $r2};
}

sub help_msg {
	print <<MSG;
Usage:
  
  perl $0 /your/standard/fastq/file/dir > list
  
  The specified directory must contain /core and /view-by-pid as two nearest
children directories.
  
  The script assumes your directory structure as:
    
    view-by-pid/\$PID/\$type/paired/\$lane/sequence/\$fastqfile
  
  Currently we only use the first \$type directory in view-by-pid/\$PID/
(e.g. if you have tumour/ and normal/, tumor/ will be ignored. This issue
will be fixed in the future.)

  The symbolic link will be resolved to the absolute path.
	
MSG
}