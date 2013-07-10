package Genome::Site::TGI::PopulationGroup;

# Adaptor for GSC PopulationGroupMember

# Do NOT use this module from anything in the GSC schema,
# though the converse will work just fine.

# This module should contain only UR class definitions,
# relationships, and support methods.

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::PopulationGroup::Member {
    table_name => '(select pg_id,member_id from population_group_member) pgm',
    id_by => [
        population_group    => { is => 'Genome::Site::TGI::PopulationGroup', id_by => 'pg_id' },
        member              => { is => 'Genome::Site::TGI::Individual', id_by => 'member_id' },
    ],
    doc => 'the association between an individual and a given population group',
    data_source => 'Genome::DataSource::Dwrac',
};

1;

