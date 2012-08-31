package Genome::Model::Build::Link;
#:adukes related code in G:M:Build

use strict;
use warnings;

use Genome;
class Genome::Model::Build::Link {
    type_name => 'genome model build link',
    table_name => 'GENOME_MODEL_BUILD_LINK',
    id_by => [
        from_build_id => { is => 'NUMBER', len => 11, implied_by => 'from_build' },
        to_build_id   => { is => 'NUMBER', len => 11, implied_by => 'to_build' },
    ],
    has => [
        role       => { is => 'VARCHAR2', len => 56 },
        from_build => { is => 'Genome::Model::Build', id_by => 'from_build_id', constraint_name => 'GMBL_FB_GMB_FK' },
        to_build   => { is => 'Genome::Model::Build', id_by => 'to_build_id', constraint_name => 'GMBL_TB_GMB_FK' },
    ],
    unique_constraints => [
        { properties => [qw/from_build_id to_build_id/], sql => 'GMBL_PK' },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;
