package Genome::Model::Tools::Sx::PhredWriter;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sx::PhredWriter {
    is => 'Genome::Model::Tools::Sx::SeqWriter',
    has => [ qual_file => { is => 'Text', is_optional => 1, }, ],
};

sub write {
    my ($self, $seq) = @_;

    Carp::confess('No sequence to write in fasta format!') if not $seq;

    $self->_write_seq($self->{_file}, $seq) or return;
    if ( $self->{_qual_file} ) {
        $self->_write_qual($self->{_qual_file}, $seq) or return;
    }

    return 1;
}

sub _write_seq {
    my ($self, $fh, $seq) = @_;

    my $header = '>'.$seq->{id};
    $header .= ' '.$seq->{desc} if defined $seq->{desc};
    $fh->print($header."\n");

    my $sequence = $seq->{seq};
    return 1 if not defined $sequence or length($sequence) == 0;

    $sequence =~ s/(.{1,60})/$1\n/g; # 60 bases per line
    $fh->print($sequence);

    return 1;
}

sub _write_qual {
    my ($self, $fh, $seq) = @_;

    $fh->print('>'.$seq->{id}."\n");

    return 1 if not defined $seq->{qual} or length($seq->{qual}) == 0;

    my $qual_string = join(' ', map { ord($_) - 33 } split('', $seq->{qual}));
    $qual_string .= ' ';
    $qual_string =~ s/((\d\d?\s){1,25})/$1\n/g;
    $qual_string =~ s/ \n/\n/g;
    $fh->print($qual_string);

    return 1;
}

1;

