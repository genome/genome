package Genome::Model::SV;

use strict;
use warnings;

use Genome;
class Genome::Model::SV {
    type_name => 'genome model sv',
    table_name => 'GENOME_MODEL_SV',
    id_by => [
        variant_id => { is => 'NUMBER', len => 10 },
    ],
    has => [
        chromosome_1    => { is => 'VARCHAR2', len => 10 },
        chromosome_2    => { is => 'VARCHAR2', len => 10 },
        event_type      => { is => 'VARCHAR2', len => 20 },
        orientation_1   => { is => 'VARCHAR2', len => 10 },
        orientation_2   => { is => 'VARCHAR2', len => 10 },
        pos_predicted_1 => { is => 'NUMBER', len => 12 },
        pos_predicted_2 => { is => 'NUMBER', len => 12 },
        sv_size         => { is => 'NUMBER', len => 12 },
    ],
    unique_constraints => [
        { properties => [qw/chromosome_1 chromosome_2 event_type orientation_1 orientation_2 pos_predicted_1 pos_predicted_2 sv_size/], sql => 'MGSV_UK' },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;
