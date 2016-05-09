package Genome::Model::Tools::Sx::Split::ByNs;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sx::Split::ByNs {
    is  => 'Genome::Model::Tools::Sx::Base',
    has => {
        number_of_ns => {
            is => 'Integer',
            doc => 'Number of Ns to split sequence.',
        },
    },
};

sub help_brief { 'Split sequences by Ns' }

sub help_detail { 'Split sequences by a set number of Ns. Does not support splitting qualites.' }
    
sub execute {
    my $self = shift;

    return if not $self->_init;

    my $reader = $self->_reader;
    while ( my $seqs = $reader->read ) {
        for my $seq ( @$seqs ) {
            $self->_split_and_write_seq($seq);
        }
    }

    return 1;
}

sub _split_and_write_seq {
    my ($self, $seq) = @_;

    my $cnt = 0;
    my $n = $self->number_of_ns;
    for my $bases ( split(/n{$n,}/i, $seq->{seq}) ) {
        $self->_writer->write([
            {
                id => join('.', $seq->{id}, ++$cnt),
                seq => $bases,
            }
        ]);
    }

    return 1;
}

1;

