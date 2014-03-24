package Genome::Model::GenotypeMicroarray::Filter::ByChromosome;

use strict;
use warnings;
use Genome;

class Genome::Model::GenotypeMicroarray::Filter::ByChromosome {
    has => [
        exclude => {
            is => 'Text',
        },
    ],
    has_transient_optional => [
        exclude_list => {
            is_many => 1,
            is => 'Text',
        },
    ],
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    $self->exclude_list([split ",", $self->exclude]);
    return $self;
}

sub filter {
    my ($self, $variant) = @_;
    if (grep {$variant->{chrom}} $self->exclude_list) {
        return 1;
    }

    return;
}

1;

