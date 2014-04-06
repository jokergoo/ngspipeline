package CO::NGSPipeline::Getopt;

use strict;
use CO::Utils;
use Getopt::Long;
use File::Basename;
use Data::Dumper;

sub new {
	my $class = shift;
	$class = ref($class) ? ref($class) : $class;
	
	my $pipeline = shift;
	
	# internal key are started with '_'
	my $opt = bless {_before => "",
	                 _after  => "",
					 _index  => 0,   # to order the arguments
					 @_}, $class;
	
	$opt->reset_opt;
	
	return $opt;
}

# text before description for each argument
sub before {
	my $opt = shift;
	my $text = shift;
	
	$opt->{_before} = $text;
}

# text after description for each argument
sub after {
	my $opt = shift;
	my $text = shift;
	
	$opt->{_after} = $text;
}

sub reset_opt {
	my $opt = shift;
	
	my $ref0 = \1;  # just mean it is a reference
	my $ref1 = \1;
	my $ref2 = \1;
	my $ref3 = \1;
	my $ref4 = \1;
	my $ref5 = \1;
	my $ref6 = \1;
	my $ref7 = \1;
	my $ref8 = \1;
	my $ref9 = \1;
	
	# order of necessary arguments are defined here
	$opt->add($ref0, "help!",      "print help message and exit");
	$opt->add($ref1, "list=s",     "sample list, containing columns which are: 1. sample name. 2. fastq file for paired end 1, should be gzipped; 3. fastq file for paired end 2, should be gzipped. If this is single end sample, put this column empty; 4. if from multiple lanes, whether they come from the same libraries or not. (any type of strings to represent category of libraries, optional. Now it is only workable for WGBS pipeline.)");
	$opt->add($ref2, "dir=s",      "working dir, default is `analysis`. Under the working dir, there are list of directories named with sample names which are called job directory for each sample.");
	$opt->add($ref3, "tool=s",     "");
	$opt->add($ref4, "sample=s",   "subset of sample ids, should seperated by ',' (no blank)");
	$opt->add($ref5, "enforce!",   "enforce to re-run pipeline from the beginning no matter they were successfully finished or not.");
	$opt->add($ref6, "filesize=i", "If size of some output files (e.g. bam files, methylation calling files) are smaller than this value, then step is terminated. Default is 1M (1024*1024). Set it to zero or non-number strings to shut down file size checking.");
	$opt->add($ref7, "test!",      "testing mode, only generate shell script files but not submit to cluster.");
	$opt->add($ref8, "prefix=s",   "prefix added to job names");
	$opt->add($ref9, "email=s",    "user's email address which log should be sent to.");
	
	return $opt;
}

sub add {
	my $opt = shift;
	
	my $ref = shift;
	my $param = shift;
	my $desc = shift;
	
	my $name = $param; $name =~s/[|!=].*$//;

	if($opt->{$name}) {
		$opt->{$name}->{'ref'} = $ref;              # reference to store values
		$opt->{$name}->{'param'} = $param;          # argument definition send to Getopt::Long
		$opt->{$name}->{'desc'} = $desc if($desc);  # description used for help message
	} else {
		$opt->{_index} ++;
		
		$opt->{$name}->{'ref'} = $ref;
		$opt->{$name}->{'param'} = $param;
		$opt->{$name}->{'desc'} = $desc;
		$opt->{$name}->{'index'} = $opt->index;
	}
	
	return $opt;
}

# delete one argument
sub del {
	my $opt = shift;
	
	my $name = shift;
	
	delete($opt->{$name});
	return $opt;
}

# get the order of current argument
sub index {
	my $opt = shift;
	
	return $opt->{_index};
}

# option names start without _
sub opt_name {
	my $opt = shift;
	
	grep {! /^_/} keys %$opt;
}


# print help messages
sub help_msg {
	my $opt = shift;
	
	print $opt->{_before} if($opt->{_before});
	
	# length of option names
	my $max_name_len = 0;
	foreach my $name ($opt->opt_name) {
		$max_name_len = length($name) if ($max_name_len < length($name));
	}
	$max_name_len += 5;  # including leading two spaces, two '-' and following one space
	
	print "Parameters:\n\n";
	foreach my $name (sort {$opt->{$a}->{index} <=> $opt->{$b}->{index}} $opt->opt_name) {
		print "  --$name";
		
		my $i_col = length($name) + 4;
		my @words = split " ", $opt->{$name}->{'desc'};
		
		while(my $w = shift @words) {
			if($i_col < $max_name_len) {
				print " " x ($max_name_len - $i_col);
			}
			if($i_col + length($w) > 70) {
				print "\n";
				print " " x $max_name_len;
				print "$w ";
				$i_col = $max_name_len + length($w) + 1;
			} else {
				print "$w ";
				$i_col += length($w) + 1;
			}
		}
		print "\n\n";
		
	}
	
	print $opt->{_after} if($opt->{_after});
}

sub getopt {
	my $opt = shift;
	
	# if there is no arguments or specified with --help
	if(scalar(@ARGV) == 0 or grep {/\b(-h|--help|-help)\b/i} @ARGV) {
		$opt->help_msg();
		exit;
	}
	
	# generate arguments for Getopt::Long
	my %param;
	foreach my $name ($opt->opt_name) {
		$param{$opt->{$name}->{'param'}} = $opt->{$name}->{'ref'};
	}
	
	GetOptions(%param) or ($opt->help_msg, exit);
	
	$opt->validate;
}

# test whether file/dir exists
sub validate {
	my $opt = shift;
	
	my $list_ref = $opt->{'list'}->{'ref'};
	my $wd_ref = $opt->{'dir'}->{'ref'};
	my $request_sampleid_ref = $opt->{'sample'}->{'ref'};
	my $tool_ref = $opt->{'tool'}->{'ref'};
	my $filesize_ref = $opt->{'filesize'}->{'ref'};

	# check whether these optinos are specified
	if(!defined($$list_ref)) {
		die "--list should be specified\n";
	}
	
	
	my %subset_samples = map { $_ => 1} split ",", $$request_sampleid_ref;
	$$filesize_ref += 0;  # enforce it as a number

	# validate sample list file
	open F, $$list_ref or die "Cannot open $$list_ref\n";  # sample list file
	my $r1;
	my $r2;
	my $sample;
	my $n_sample = 0;
	while(my $line = <F>) {
		$line =~s/\s*$//;
		next if($line =~/^\s*$/);
		next if($line =~/^#/);
		
		my @tmp = split "\t", $line;
		
		$tmp[1] = to_abs_path($tmp[1]);
		$tmp[2] = to_abs_path($tmp[2]);
		
		unless(-e $tmp[1]) {
			die "[ $tmp[0] ] cannot find $tmp[1]\n";
		}
		unless(-e $tmp[2]) {
			die "[ $tmp[0] ] cannot find $tmp[2]\n";
		}

		
		if(basename($tmp[1]) eq basename($tmp[2])) {
			die "[ $tmp[0] ] two fastq files have the same name! check your file!\n";
		}
		
		# if specified with --sample
		if(scalar(%subset_samples) and !$subset_samples{$tmp[0]}) {
			print "$tmp[0] is not in --sample, skip this sample.\n";
			next;
		}
		
		# if no record for this sample, initialize the array reference
		if(! defined($sample->{$tmp[0]})) {
			$sample->{$tmp[0]} = {};
			$sample->{$tmp[0]}->{r1} = [];
			$sample->{$tmp[0]}->{r2} = [];
			$sample->{$tmp[0]}->{library} = [];
		}
		
		push(@{$sample->{$tmp[0]}->{r1}}, $tmp[1]);
		push(@{$sample->{$tmp[0]}->{r2}}, $tmp[2]);
		
		push(@{$sample->{$tmp[0]}->{library}}, defined($tmp[3]) ? $tmp[3] : "default");
		
		$n_sample ++;
	}

	# calculate total file size for all fastq files
	foreach my $sid (keys %$sample) {
		$sample->{$sid}->{sum_filesize} = 0;
		foreach my $f (@{$sample->{$sid}->{r1}}) {
			$sample->{$sid}->{sum_filesize} += -s $f;
		}
		foreach my $f (@{$sample->{$sid}->{r2}}) {
			$sample->{$sid}->{sum_filesize} += -s $f;
		}
	}

	if(!defined($$wd_ref)) {
		$$wd_ref = "analysis";
	}
	$$wd_ref = to_abs_path($$wd_ref);
	# seems set mode of the dir to 0755 not always successful
	-e $$wd_ref ? 1: mkdir $$wd_ref, 0775 || die "cannto create dir: $$wd_ref with mode 0775\n";

	$$wd_ref =~s/\/$//;
	
	$$tool_ref = lc($$tool_ref);
	
	print "Working directory is $$wd_ref.\n";
	print "Totally $n_sample samples with ". scalar(keys %$sample)." unique sample ids \n";
	print "Using $$tool_ref.\n\n";

	$$list_ref = $sample;
}

1;
