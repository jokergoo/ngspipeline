package CO::Utils;

###########################################
# useful util subroutines
###########################################

use strict;
use File::Spec;
use Data::Dumper;
use File::Basename;
use Cwd 'abs_path';

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(to_abs_path
                 dir_tree_as_hash);

sub to_abs_path {
	my $path = shift;
	return File::Spec->file_name_is_absolute($path) ? $path : File::Spec->rel2abs($path);
}


sub dir_tree_as_hash {
	my $dir = shift;
	$dir = to_abs_path($dir);
	if(-f $dir) {
		return {$dir => 1}
	} else {
		return children_dir($dir);
	}
}

sub children_dir {
	my $parent = shift;
	
	my $hash;
	my @children = glob("$parent/*");
	foreach my $child (@children) {
		
		if(-d $child) {
			$hash->{basename($child)} = children_dir($child);
		} elsif(-l $child) {
			$hash->{$child} = 1
		} else {
			$hash->{$child} = 1;
		}
	}
	return $hash;
}

sub replace_suffix {
	my $filename = shift;
	my $old = shift;
	my $new = shift;
	
	$filename =~s/$old$/$new/;
	
	return $filename;
}

1;
