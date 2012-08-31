package Genome::Utility::FastaStreamOut;

#:eclark 11/17/2009 Code review.

# One of many reader/writer classes in Genome/Utility.  Need a consistant naming convention and interface for all of them.

use strict;
use warnings;
use Data::Dumper;

sub new{
    my $class = shift;
    my $io = shift;
    my $linelength = shift;
    $linelength||=60;
    die "Need IO::Handle" unless $io->isa('IO::Handle');
    my $self = bless({_io => $io, current_line_avail => $linelength, linelength=>$linelength},$class);
    return $self;
}

sub print_header{
    my ($self, $header) = @_;
    #don't print leading newline if we're at the top of the file
    $self->{_io}->print("\n") unless $self->{current_line_avail} == $self->{linelength};
    $self->{_io}->print("$header\n") or $self->fatal_msg("can't write header $header");
    $self->{current_line_avail} = $self->{linelength};
    return 1;
}

sub print{
    my $self = shift;
    my $avail = $self->{current_line_avail};
    my $io = $self->{_io};
    while ($_ = shift @_) {
        next unless $_;
        my $next = substr($_,0,$avail);
        $io->print($next);
        $avail -= length($next);
        if ($avail == 0) {
            $io->print("\n");
            $avail = $self->{linelength};
        }                    
        $_ = substr($_,length($next));
        redo if length($_);        
    }
    $self->{current_line_avail}=$avail;
    return 1;
}

sub close{
    my $self = shift;
    my $io = $self->{_io};
    $io->print("\n");
    $io->close();
}

=pod

=head1 FastaStreamOut
Simple output writer for taking arbitrary length input and writing output of max line length(60)

my $ob = Genome::Utility::FastaStreamOut->new(<file>);

$ob->print_header(">Sequence_1");

while{
    ...
    (create $seq of arbitrary length
    ...
    $ob->print("$seq);
}

=head2 Subs

=head3 print_header($string)
prints $string followed by a newline to the file

=head3 print($string)
prints $string to the file, automatically inserting a newline when the max line_length(60) is reached

=cut

1;
