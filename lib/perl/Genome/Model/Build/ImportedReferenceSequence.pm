package Genome::Model::Build::ImportedReferenceSequence;

use Genome;

use strict;
use warnings;


class Genome::Model::Build::ImportedReferenceSequence {
    is => 'Genome::Model::Build::ReferenceSequence',
    has => [
        species_name => { via => 'subject', to => 'name' },

        derived_from => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            doc => 'Identifies the parent build from which this one is derived, if any.',
        },
        coordinates_from => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            doc => 'Used to indicate that this build is on the same coordinate system as another.',
        },
        append_to => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            doc => 'If specified, the created reference will be logically appended to the one specified by this parameter for aligners that support it.',
        },
        combines => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            doc => 'If specified, merges several other references into one.', 
        },
    ],
};


sub get{
    my $self = shift;
    my @results = $self->SUPER::get(@_);
    return $self->SUPER::get(@_) if @results;

    return Genome::Model::Build::ReferenceSequence->get(@_);
}


1;
