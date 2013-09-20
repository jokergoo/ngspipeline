package CO::FastQ;

use CO::FastQ::Read;
use strict;

=pod

=head1 NAME

CO::FastQ - Simple FastQ class

=head1 SYNOPSIS

  use CO::FastQ;
  
  my $fastq;
  $fastq = CO::FastQ->new(file => "r1.fastq");
  $fastq = CO::FastQ->new(file => "r1.fasq.gz");
  
  open my $fh, ">", "r1.fastq";
  $fastq = CO::FastQ->new(file => $fh);
  open my $fh, "zcat -c r1.fastq.gz | ";
  $fastq = CO::FastQ->new(file => $fh);
  
  while(my $read = $fastq->next) {
	print "read #".$fastq->i."\n";
  }

=head1 DESCRIPTION

This is a simple class for FastQ files. FastQ file contains short sequence records
repeatly with block of 4 lines. The class in fact stores file handle of the FastQ
file and the index of read that has been read.

=head2 Subroutines
  
=cut


=pod

=over 4

=item C<CO::FastQ-E<gt>new(file =E<gt> $file)>

C<$file> can either be uncompressed FastQ file, gzipped FastQ file (interpreted
from file name) or an already opened file handle. The function will transform
files into file handles.

=cut
sub new {
	my $class = shift;
	$class = ref($class) ? ref($class) : $class;
	
	my $self = {"file" => undef,
	            "i"    => 0,
	            @_};
	
	if(!defined($self->{file})) {
		die "`file` should be specified.\n"
	}
	
	if(ref($self->{file}) eq "GLOB") {  # already opened file handle
	
		$self->{fh} = $self->{file};
		#print STDERR "file($self->{file}) is a file handle\n";
		
	} elsif($self->{file} =~/\.gz$/) {  # gzipped fastq file
	
		open my $fh, "zcat -c $self->{file} | " or die "cannot open pipe `zcat -c $self->{file} | `\n";
		$self->{fh} = $fh;
		#print STDERR "file($self->{file}) is gzipped file\n";
		
	} else {  # uncompressed fastq file
	
		open my $fh, $self->{file} or die "cannot open $self->{file}\n";
		$self->{fh} = $fh;
		#print STDERR "file($self->{file}) is uncompressed fastq file\n";
		
	} 
	
	return bless $self, $class;
}

=pod

=item C<$fastq-E<gt>next>

return next read. It returns a L<CO::FastQ::Read> object. If it reaches the end
of FastQ file, it returns C<undef>.

=cut
sub next {
	my $self = shift;
	
	my $fh = $self->{fh};
	
	# EOF
	if(eof($fh)) {
		return undef;
	}
	
	my $line1 = <$fh>; chomp $line1;
	my $line2 = <$fh>; chomp $line2;
	my $line3 = <$fh>; chomp $line3;
	my $line4 = <$fh>; chomp $line4;
	
	# for efficiency, we store the reference.
	my $read = CO::FastQ::Read->new(line1 => \$line1,
                                    line2 => \$line2,
                                    line3 => \$line3,
                                    line4 => \$line4,
                                    start => 0,
                                    end => length($line2),);
	
	$self->{i} ++;
	
	return $read;
}

=pod

=item C<$fastq-E<gt>i>

Index of read that has beem read.

=cut
sub i {
	my $self = shift;
	$self->{i};
}

=pod

=back

=head1 AUTHOR

Zuguang Gu E<lt>z.gu@dkfz.deE<gt>

=cut

1;
