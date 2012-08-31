package Genome::Model::Link;
#:adukes related code in G:Model

use strict;
use warnings;

use Genome;
class Genome::Model::Link {
    type_name => 'genome model link',
    table_name => 'GENOME_MODEL_LINK',
    id_by => [
        from_model_id => { is => 'NUMBER', len => 11, implied_by => 'from_model' },
        to_model_id   => { is => 'NUMBER', len => 11, implied_by => 'to_model' },
    ],
    has => [
        role       => { is => 'VARCHAR2', len => 56 },
        from_model => { is => 'Genome::Model', id_by => 'from_model_id', constraint_name => 'GML_FB_GM_FK' },
        to_model   => { is => 'Genome::Model', id_by => 'to_model_id', constraint_name => 'GML_TB_GM_FK' },
    ],
    unique_constraints => [
        { properties => [qw/from_model_id to_model_id/], sql => 'GML_PK' },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;
