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
	die "'view-by-pid/' should be the first level child directory under $ARGV[0].\n";
}

my $match = $ARGV[1] || '.*';
$match = qr($match);

my $list = read_std_dir($ARGV[0], $match);

my $pid = $list->{pid};
my $r1 = $list->{r1};
my $r2 = $list->{r2};



for(my $i = 0; $i < len($pid); $i ++) {
	print "$pid->[$i]\t$r1->[$i]\t$r2->[$i]\n";
}

sub read_std_dir {
	my $dir = shift;
	my $match = shift;
	
	$dir =~s/\/$//;
	
	if(! -d $dir) {
		help_msg();
		die "cannot find dir: $dir";
	} elsif(!($dir=~/view-by-pid$/)) {
		help_msg();
		die "dir should end with view-by-pid\n";
	}
	
	my $tree = dir_tree_as_hash($dir);
	
	my $pid_dir = [ sort grep {-d "$dir/$_"} keys %{$tree} ];
	
	my $r1 = [];
	my $r2 = [];
	my $pid = [];

	for(my $i = 0; $i < len($pid_dir); $i ++) {
		
		my $type_dir = [ sort grep {-d "$dir/$pid_dir->[$i]/$_"} keys %{$tree->{$pid_dir->[$i]}} ];
		
		for(my $j = 0; $j < len($type_dir); $j ++) {
			my @lanes = sort grep {-d "$dir/$pid_dir->[$i]/$type_dir->[$j]/paired/$_"} keys %{$tree->{$pid_dir->[$i]}->{$type_dir->[$j]}->{paired}};
			my @lanes2 = sort grep {-d "$dir/$pid_dir->[$i]/$type_dir->[$j]/single/$_"} keys %{$tree->{$pid_dir->[$i]}->{$type_dir->[$j]}->{single}};
			push(@lanes, @lanes2);
			
			for(my $k = 0; $k < scalar(@lanes); $k ++) {
				if($lanes[$k] =~/^run/) {
					my @reads =  sort grep {/$match/} keys %{$tree->{$pid_dir->[$i]}->{$type_dir->[$j]}->{paired}->{$lanes[$k]}->{sequence}};
					if(scalar(@reads)) {
						for(my $x = 0; $x < scalar(@reads); $x +=2) {
							push(@$pid, "$pid_dir->[$i]_$type_dir->[$j]");
							push(@$r1, $reads[$x]);
							push(@$r2, $reads[$x+1]);
						}
					}
					
					my @reads =  sort grep {/$match/} keys %{$tree->{$pid_dir->[$i]}->{$type_dir->[$j]}->{single}->{$lanes[$k]}->{sequence}};
					if(scalar(@reads)) {
						for(my $x = 0; $x < scalar(@reads); $x ++) {
							push(@$pid, "$pid_dir->[$i]_$type_dir->[$j]");
							push(@$r1, $reads[$x]);
							push(@$r2, "");
						}
					}
				}
			}
		}
	}
	return {pid => $pid, r1 => $r1, r2 => $r2};
}

sub help_msg {
	print <<MSG;
Generate sample list file from directory with standard structure
	
Usage:
  
  perl $0 /path/to/view-by-pid > list
  
  The specified directory must end with 'view-by-pid'
  
  The script assumes your directory structure as:
    
    view-by-pid/\$PID/\$type/paired|single/\$lane/sequence/\$r1_fastq
    view-by-pid/\$PID/\$type/paired|single/\$lane/sequence/\$r2_fastq

  The sample name would be \$PID_\$type
  
  The second option can be regular expression to match the final
FastQ files. By default we assume there are only two files under
sequence/ folder. But if there is more files, you can specify
this option to match file names.

  E.g. if there are two files r1.fq.gz and r2.fq,gz, you can write
the command as:

  perl $0 /path/to/view-by-pid 'r[12]\\.fq\\.gz\$' > list
  
  DO NOT write slash around the regular expression.
  
MSG
}