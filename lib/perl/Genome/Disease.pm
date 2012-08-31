package Genome::Disease;

use strict;
use warnings;

use Genome;
class Genome::Disease {
    type_name => 'disease',
    table_name => 'DISEASE',
    id_by => [
        id => { is => 'NUMBER', len => 20 },
    ],
    has => [
        common_name            => { is => 'VARCHAR2', len => 20, is_optional => 1 },
        disease_parent_disease => { is => 'Genome::Disease', id_by => 'parent_disease_id', constraint_name => 'DISEASE_PARENT_FK', is_optional => 1 },
        full_name              => { is => 'VARCHAR2', len => 256 },
        parent_disease_id      => { is => 'NUMBER', len => 20, implied_by => 'disease_parent_disease', is_optional => 1 },
        population_group_id    => { is => 'NUMBER', len => 15, is_optional => 1 },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    doc => 'One row per disease which is an explicit subject of study.',
};

1;
