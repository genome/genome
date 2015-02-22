package Genome::Model::Build::Link;
#:adukes related code in G:M:Build

use strict;
use warnings;

use Genome;
class Genome::Model::Build::Link {
    table_name => 'model.build_link',
    type_name => 'genome model build link',
    id_by => [
        from_build_id => { is => 'Text', len => 64 },
        to_build_id => { is => 'Text', len => 64 },
    ],
    has => [
        role => { is => 'Text', len => 56 },
        from_build => {
            is => 'Genome::Model::Build',
            id_by => 'from_build_id',
            constraint_name => 'GMBL_FB_GMB_FK',
        },
        to_build => {
            is => 'Genome::Model::Build',
            id_by => 'to_build_id',
            constraint_name => 'GMBL_TB_GMB_FK',
        },
    ],
    unique_constraints => [
        { properties => [qw/from_build_id to_build_id/], sql => 'GMBL_PK' },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;
