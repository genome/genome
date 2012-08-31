package Genome::Model::Build::ImportedReferenceSequence;
use strict;
use warnings;
use Genome;

# all of the logic has been moved into a new class
# once fully stable, we will do a db update to flip the class names in the model/pp/build tables.

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
