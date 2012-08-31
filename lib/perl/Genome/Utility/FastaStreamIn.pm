package Genome::Utility::FastaStreamIn;

#:eclark 11/17/2009 Code review.

# One of many reader/writer classes in Genome/Utility.  Need a consistant naming convention and interface for all of them.

use strict;
use warnings;
use Data::Dumper;

sub new{
    my $class = shift;
    my $io = shift;
    die "Need IO::Handle" unless $io->isa('IO::Handle');
    my $next_line = $io->getline;
    my $self = bless({_io => $io, next_line => $next_line },$class);
    return $self;
}

sub next_line{
    my $self = shift;
    return unless $self->{next_line}; 
    if (substr($self->{next_line},0,1) eq '>'){
        return undef;
    }
    my $next_line = $self->{next_line};
    $self->{next_line} = $self->{_io}->getline;
    chomp $next_line;
    return uc $next_line;
}

sub next_header{
    my $self = shift;
    if (not defined $self->{next_line}){
        return undef;
    }elsif (!substr($self->{next_line},0,1) eq '>' ){
        while($self->next_line){}
        return unless defined $self->{next_line};
    }
    $self->{current_header} = $self->{next_line};
    $self->{next_line} = $self->{_io}->getline;
    return $self->{current_header};
}

sub current_header{
    my $self = shift;
    my $line = $self->{current_header};
    chomp $line;
    return $line;
}

sub _lookahead{ #TODO rewrite this to handle all types of streams, not just files
    my ($self, $distance) = @_;
    my $pos = $self->{_io}->tell;
    my $string = $self->{next_line};
    chomp $string;
    while (length $string < $distance){
        my $next_line =$self->{_io}->getline;
        last if substr($next_line,0,1) eq '>'; 
        chomp $next_line;
        $string.=$next_line;
    }
    $self->{_io}->seek($pos, 0);
    return substr($string, 0, $distance);
}

=pod

=head1 FastaStreamIn
fasta file input stream used in genome-model-tools apply-diff-to-fasta

=head2 Synopsis

This streams through a fasta file, returning fasta sequence and header data;

my $fs = Genome::Utility::FastaStreamIn->new( <file_name> );
my $first_header_line = $fh->next_header;
my $sequence;
while (my $line = $fh->next_line){
$sequence.=$line;
}

This object can be used to read a fasta file in a streaming format. If you want to recreate the header line, use current_header() before advancing to the next fasta section

=head2 Subs

=head3 next_header
reads and returns the next header line in the file.  If the current fasta sequence section hasn't been fully streamed through, the file is adavanced until the end of the current fasta section.

=head3 next_line
reads and returns the next line in the current fasta section.  Returns undef when the end of the section is reached. 

=head3 current_header
Returns the fasta header line of the current section. Does not advance through the file.

=head3 lookahead(int $distance)
returns the next $distance chars of the current fasta section, from the current position in that section.  If the lookahead distance exceeds the length of the fasta section, returns the sequence until the end of the section.   Does not advance the file.
=cut

1;
