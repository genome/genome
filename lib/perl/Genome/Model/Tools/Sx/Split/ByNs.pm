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

    my $init = $self->_init;
    return if not $init;

    my $reader = $self->_reader;
    my $writer = $self->_writer;
    while ( my $seqs = $reader->read ) {
        for my $seq ( @$seqs ) {
            $writer->write( $self->_split_seq($seq) );
        }
    }

    return 1;
}

sub _split_seq {
    my ($self, $seq) = @_;

    my @seqs;
    my $cnt = 0;
    my $n = $self->number_of_ns;
    for my $bases ( split(/n{$n,}/i, $seq->{seq}) ) {
        push @seqs, {
            id => join('.', $seq->{id}, ++$cnt),
            seq => $bases,
        };
    }

    return \@seqs;
}

1;

