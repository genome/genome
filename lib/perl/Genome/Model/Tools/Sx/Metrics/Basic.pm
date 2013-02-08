package Genome::Model::Tools::Sx::Metrics::Basic;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sx::Metrics::Basic {
    is => 'Genome::Model::Tools::Sx::Metrics',
    has => [
        (map { 
            $_ => {
                #is => 'UR::Value',
                default_value => 0,
            }
        } __PACKAGE__->metric_names()),
    ],
};


sub metric_names { return (qw/ bases count /); }

sub add_sequence {
    my ($self, $seq) = @_;

    $self->bases( $self->bases + length($seq->{seq}) );
    $self->count( $self->count + 1 );

    return 1;
}

sub add_sequences {
    my ($self, $seqs) = @_;

    for my $seq ( @$seqs ) {
        $self->add_sequence($seq);
    }

    return 1;
}

1;

