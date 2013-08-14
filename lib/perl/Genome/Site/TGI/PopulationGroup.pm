package Genome::Site::TGI::PopulationGroup;

# Adaptor for GSC PopulationGroup

# Do NOT use this module from anything in the GSC schema,
# though the converse will work just fine.

# This module should contain only UR class definitions,
# relationships, and support methods.

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::PopulationGroup {
    is => 'Genome::Site::TGI::Measurable',
    table_name => 'POPULATION_GROUP',
    id_by => [
        individual_id => { is => 'Number', len => 10, column_name => 'PG_ID' },
    ],
    has => [
        name => { is => 'Text', len => 64, doc => 'Name of the Population Group.', },
        taxon => { is => 'Genome::Site::TGI::Taxon', id_by => 'taxon_id', doc => 'The taxon to which this individual belongs', },
        species_name => { via => 'taxon' },
        description => { is => 'Text', is_optional => 1, len => 500, doc => 'Description', },
        subject_type => { is => 'Text', is_constant => 1, value => 'population group', column_name => '', },
    ],
    has_many => [
        member_links        => { is => 'Genome::Site::TGI::PopulationGroup::Member', reverse_id_by => 'population_group' },
        members             => { is => 'Genome::Site::TGI::Individual', via => 'member_links', to => 'member' },
        samples => { 
            is => 'Genome::Site::TGI::Sample', 
            is_many => 1,
            reverse_id_by => 'source',
        },
        sample_names => {
            via => 'samples',
            to => 'name', is_many => 1,
        },
    ],
    doc => 'an defined, possibly arbitrary, group of individual organisms',
    data_source => 'Genome::DataSource::Dwrac',
};

sub common_name { # not in the table, but exepected by views
    return $_[0]->name;
}

1;

