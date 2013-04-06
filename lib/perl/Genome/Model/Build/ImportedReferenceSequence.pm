package Genome::Model::Build::ImportedReferenceSequence;
use strict;
use warnings;
use Genome;

# all of the logic has been moved into a new class
# once fully stable, we will do a db update to flip the class names in the model/pp/build tables.

sub get{
    my $self = shift;
    my @results = $self->SUPER::get(@_);
    return $self->SUPER::get(@_) if @results;

    return Genome::Model::Build::ReferenceSequence->get(@_);
}

# the class def is here, to avoid circularity issues
require Genome::Model::ReferenceSequence;

1;

