package Genome::Model::Build::ImportedReferenceSequence;

use Genome;

use strict;
use warnings;


class Genome::Model::Build::ImportedReferenceSequence {
    is => 'Genome::Model::Build::ReferenceSequence',
};


sub get{
    my $self = shift;
    my @results = $self->SUPER::get(@_);
    return $self->SUPER::get(@_) if @results;

    return Genome::Model::Build::ReferenceSequence->get(@_);
}

1;
