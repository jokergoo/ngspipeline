#!/bin/perl/bin/

BEGIN {
	use File::Basename;
	unshift(@INC, dirname($0)."/lib");
}

use strict;
use CO::Utils;
use File::Basename;
use List::Vectorize;

use Data::Dumper;

if(! @ARGV) {
	help_msg();
	exit 0;
}

unless( grep {~/^view-by-pid$/} split "\n", `ls $ARGV[0]` ){
	die "'view-by-pid' should be the first level child directory under $ARGV[0].\n";
}

my $list = read_std_dir($ARGV[0]);
my $pid = $list->{pid};
my $type = $list->{tpye};
my $r1 = $list->{r1};
my $r2 = $list->{r2};

my $need_type = 0;
if(any(sapply($type, sub { len($_[0]) > 1 }))) {
	$need_type = 1;
}

for(my $i = 0; $i < len($pid); $i ++) {
	for(my $j = 0; $j < len($tpye->[$i])) {
		for(my $k = 0; $k < len($r1->[$i]->[$j]); $j ++) {
			if($need_type) {
				print "$r1->[$i]->[$j]->[$k]\t$r2->[$i]->[$j]->[$k]\t$pid->[$i]_$type->[$i]->[$j]\n";
			} else {
				print "$r1->[$i]->[$j]->[$k]\t$r2->[$i]->[$j]->[$k]\t$pid->[$i]\n";
			}
		}
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
		die "cannot find view-by-pid/ directory in $dir\n";
	}
	
	my $pid = [ keys %{$tree->{"view-by-pid"}} ];
	my $r1 = [];
	my $r2 = [];
	my $type = [];
	
	for(my $i = 0; $i < len($pid); $i ++) {
		my $type->[$i] = [ keys %{$tree->{"view-by-pid"}->{$pid->[$i]}} ];
		
		for(my $j = 0; $j < len($tpye->[$i]); $j ++) {
			my @lanes = keys %{$tree->{"view-by-pid"}->{$pid->[$i]}->{$type->[$i]->[$j]}->{paired}};
			
			for(my $k = 0; $k < scalar(@lanes); $k ++) {
				($r1->[$i]->[$j]->[$k], $r2->[$i]->[$j]->[$k]) = keys %{$tree->{"view-by-pid"}->{$pid->[$i]}->{$type->[$i]->[$j]}->{paired}->{$lanes[$k]}->{sequence}};
			}
		}
	}
	return {pid => $pid, tpye => $type, r1 => $r1, r2 => $r2};
}

sub help_msg {
	print <<MSG;
Usage:
  
  perl $0 /your/standard/fastq/file/dir > list
  
  The specified directory must contain core/ and view-by-pid/ as two nearest
child directories.
  
  The script assumes your directory structure as:
    
    view-by-pid/\$PID/\$type/paired/\$lane/sequence/$r1_fastq
	view-by-pid/\$PID/\$type/paired/\$lane/sequence/$r2_fastq
	
MSG
}