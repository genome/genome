package Genome::Site::TGI::Individual;

# Adaptor for GSC Organism Individual

# Do NOT use this module from anything in the GSC schema,
# though the converse will work just fine.

# This module should contain only UR class definitions,
# relationships, and support methods.

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::Individual {
    is => 'UR::Object',
    table_name => 'ORGANISM_INDIVIDUAL',
    id_by => [
        individual_id => { is => 'Number', len => 10, column_name => 'ORGANISM_ID' },
    ],
    has => [
        name => { is => 'Text', len => 64, column_name => 'FULL_NAME', doc => 'Name of the individual', },
        taxon => { is => 'Genome::Site::TGI::Taxon', id_by => 'taxon_id', doc => 'The taxon to which this individual belongs', },
        species_name    => { is => 'Text', via => 'taxon' },
        description => { is => 'Text', is_optional => 1, len => 500, doc => 'Description', },
        subject_type => { is => 'Text', is_constant => 1, value => 'organism individual', column_name => '', },
    ],
    has_optional => [
        father  => { is => 'Genome::Site::TGI::Individual', id_by => 'father_id', doc => 'Father of this individual', },
        father_name => { is => 'Text', via => 'father', to => 'name' },
        mother  => { is => 'Genome::Site::TGI::Individual', id_by => 'mother_id', doc => 'Mother of this individual', },
        mother_name => { is => 'Text', via => 'mother', to => 'name' },
        upn => { 
            is => 'Text', 
            column_name => 'NAME',
            doc => 'fully qualified internal name for the patient', 
        },
        common_name => { 
            is => 'Text',
            len => 10,
            doc => 'a name like "aml1" for the patient, by which the patient is commonly referred-to in the lab' 
        },
        gender => { 
            is => 'Text',
            len => 16,
            doc => 'when the gender of the individual is known, this value is set to male/female/...' 
        },
        ethnicity       => { 
            is => 'Text',
            len => 64,
            doc => 'the "ethnicity" of the individual, Hispanic/Non-Hispanic/...'
        },
        race => { 
            is => 'Text',
            len => 64,
            doc => 'the "race" of the individual, African American/Caucasian/...'
        },
        nomenclature => {
            is => 'Text',
            len => 64,
            default_value => 'WUGC',
            doc => 'Nomenclature',
        },
        samples => { 
            is => 'Genome::Site::TGI::Sample', 
            is_many => 1,
            reverse_id_by => 'source',
        },
        sample_names => {
            is => 'Text',
            via => 'samples',
            to => 'name',
            is_many => 1,
        },
    ],
    data_source => 'Genome::DataSource::Dwrac',
};

sub __display_name__ {
    return $_[0]->name.' ('.$_[0]->id.')';
}

1;

