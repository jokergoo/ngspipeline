package CO::FastQ::Read;


use strict;
use base 'CO::FastQ::Trim';



=pod

=head1 NAME

CO::FastQ::Read - Simple Read class

=head1 SYNOPSIS

  use CO::FastQ::Read;
  
  my $read = CO::FastQ::Read->new{line1 => \"\@read1",
                                  line2 => \"CGAGTCGTACGTAGTCGTAC",
                                  line3 => \"+",
                                  line4 => \"IIIIIIIIIIIIIIIIIIII",
                                  begin => 0,
                                  end => length("CGAGTCGTACGTAGTCGTAC")};
                                  
  print $read->header;
  print $read->name;
  print $read->seq;
  print $read->qual_str;
  print $read->record;
  print $read->length;
  print join ",", @{$read->qual};

=head1 DESCRIPTION

This is a simple class for FastQ read. For efficiency, four lines in a read are
stored as a scalar reference. Also the start and end position of the sequence
are stored. If trimming is applied, start position and end position are modified
while the origin seqeuence are retained.

=head2 Subroutines
  
=cut


=pod

=over 4

=item C<CO::FastQ::Read-E<gt>new(HASH)>

Initialize the C<CO::FastQ::Read> object. Normally it is created in the
C<CO::FastQ> object.

=cut
sub new {
	my $class = shift;
	$class = ref($class) ? ref($class) : $class;
	
	my $self = {"line1" => undef,
	            "line2" => undef,
		        "line3" => undef,
		        "line4" => undef,
		        "begin" => undef,
 		        "end" => undef,
	            @_};
	
	return bless $self, $class;
}

=pod

=item C<$read-E<gt>header>

The first line in the read record

=cut
sub header {
	my $self = shift;
	${$self->{line1}};
}

=pod

=item C<$read-E<gt>name>

Read name extracted from the first line in the read record

=cut
sub name {
	my $self = shift;
	my $header = $self->header;
	my $i;
	for($i = 0; $i < length($header); $i ++) {
		if(substr($header, $i, 1) eq " ") {
			last;
		}
	}
	return substr($header, 0, $i);
}

=pod

=item C<$read-E<gt>seq>

Short sequence substringed by the start position and end position

=cut
sub seq {
	my $self = shift;
	substr(${$self->{line2}}, $self->{begin}, $self->{end} - $self->{begin});
}

=pod

=item C<$read-E<gt>qual_str>

quality string substringed by the start position and end position

=cut
sub qual_str {
	my $self = shift;
	substr(${$self->{line4}}, $self->{begin}, $self->{end} - $self->{begin});
}

=pod

=item C<$read-E<gt>qual>

numeric quality values

=cut
sub qual {
	my $self = shift;
	my $base = shift || 33;
	
	my @let = split "", $self->qual_str;
	[ map {ord($_) - $base} @let ];

}

=pod

=item C<$read-E<gt>record>

Totally four lines

=cut
sub record {
	my $self = shift;
	
	return "${$self->{line1}}\n".$self->seq."\n${$self->{line3}}\n".$self->qual_str."\n";
	
	
}

=pod

=item C<$read-E<gt>length>

Length of the sequence

=cut
sub length {
	my $self = shift;
	
	return $self->{end} - $self->{begin};
}


=pod

=back

=head1 AUTHOR

Zuguang Gu E<lt>z.gu@dkfz.deE<gt>

=cut
1;
