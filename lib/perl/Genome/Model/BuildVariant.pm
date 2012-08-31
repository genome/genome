# Review: gsanders - We probably need to add somatic_score and mapping_score to this table. 
# It nice to have a simple bridge table but we need these build specific stats at some point unless there is a better way

package Genome::Model::BuildVariant;

use strict;
use warnings;

use Genome;
class Genome::Model::BuildVariant {
    type_name => 'genome model build variant',
    table_name => 'GENOME_MODEL_BUILD_VARIANT',
    er_role => 'bridge',
    id_by => [
        build_id   => { is => 'NUMBER', len => 10, implied_by => 'build' },
        variant_id => { is => 'NUMBER', len => 10, implied_by => 'variant' },
    ],
    has => [
        build           => { is => 'Genome::Model::Build', id_by => 'build_id', constraint_name => 'GMBV_GMB_FK' },
        variant         => { is => 'Genome::Model::Variant', id_by => 'variant_id', constraint_name => 'GMBV_GMV_FK' },
        mapping_quality => { is => 'NUMBER', len => 5, is_optional => 1 },
        somatic_quality => { is => 'NUMBER', len => 5, is_optional => 1 },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;
