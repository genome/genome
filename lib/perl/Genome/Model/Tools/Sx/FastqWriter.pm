package Genome::Model::Tools::Sx::FastqWriter;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sx::FastqWriter {
    is => 'Genome::Model::Tools::Sx::SeqWriter',
};

sub write {
    my ($self, $seq) = @_;

    Carp::confess('No sequence to write in fastq format!') if not $seq;
    Carp::confess('Sequence and qulity lengths are not the same!'.Data::Dumper::Dumper($seq)) if length($seq->{seq}) != length($seq->{qual});

    $self->{_file}->print(
        join(
            "\n",
            '@'.$seq->{id}.( defined $seq->{desc} && $seq->{desc} ne '' ? ' '.$seq->{desc} : '' ),
            $seq->{seq},
            '+',
            $seq->{qual},
        )."\n"
    );

    return 1;
}

1;

