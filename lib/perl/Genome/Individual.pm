package Genome::Individual;

use strict;
use warnings;

use Genome;
use Carp;

class Genome::Individual {
    is => ['Genome::SampleSource','Genome::Searchable'],
    has => [
        individual_id => {
            calculate_from => 'id',
            calculate => q{ return $id; },
        },
        subject_type => {
            is_constant => 1,
            is_classwide => 1,
            value => 'sample_group',
        },
        taxon_id => {
            is => 'Number',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'taxon_id' ],
            is_mutable => 1,
            doc => 'Taxon ID for this individual',
        },
        taxon => {
            is => 'Genome::Taxon',
            id_by => 'taxon_id',
        },
        species_name => {
            via => 'taxon'
        },
    ],
    has_optional => [
        father_id => {
            is => 'Number',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'father_id' ],
            is_mutable => 1,
            doc => 'ID of this individual\' father',
        },
        father => {
            is => 'Genome::Individual',
            id_by => 'father_id',
            doc => 'Genome::Individual object that represents this individual\'s father',
        },
        father_name => {
            via => 'father',
            to => 'name',
            doc => 'Name of this individual\'s father',
        },
        mother_id => {
            is => 'Number',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'mother_id' ],
            is_mutable => 1,
            doc => 'ID of this individual\'s mother',
        },
        mother => {
            is => 'Genome::Individual',
            id_by => 'mother_id',
            doc => 'Your mom',
        },
        mother_name => {
            via => 'mother',
            to => 'name',
            doc => 'Name of this individual\'s mother',
        },
        upn => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'upn' ],
            is_mutable => 1,
            doc => 'Fully qualified internal name for the patient',
        },
        common_name => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'common_name' ],
            is_mutable => 1,
            doc => 'A name like "aml1" for the patient, by which the patient is commonly referred-to in the lab'
        },
        gender => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'gender' ],
            is_mutable => 1,
            doc => 'when the gender of the individual is known, this value is set to male/female/...'
        },
        ethnicity => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'ethnicity' ],
            is_mutable => 1,
            doc => 'The "ethnicity" of the individual, Hispanic/Non-Hispanic/...'
        },
        race => {
            is => 'Text',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'race' ],
            is_mutable => 1,
            doc => 'The "race" of the individual, African American/Caucasian/...'
        },
    ],
    has_many_optional => [
        samples => {
            is => 'Genome::Sample',
            reverse_id_by => 'source',
            doc => 'Sample extracted from this individual',
        },
        sample_names => {
            is => 'Text',
            via => 'samples',
            to => 'name',
        },
    ],
};

sub get_source {
    my $self = shift;
    return $self->taxon;
}

1;

