# Review: gsanders - It sucks having extra columns in the bridge table but there seems no better way to store build specific information in a structured way
# This deserves some thought. The only reason BuildVariant and BuildSV are separate is this build specific information.

package Genome::Model::BuildSV;

use strict;
use warnings;

use Genome;
class Genome::Model::BuildSV {
    type_name => 'genome model build sv',
    table_name => 'GENOME_MODEL_BUILD_SV',
    id_by => [
        build_id   => { is => 'NUMBER', len => 10, implied_by => 'build' },
        variant_id => { is => 'NUMBER', len => 10, implied_by => 'sv' },
    ],
    has => [
        sv                => { is => 'Genome::Model::SV', id_by => 'variant_id', constraint_name => 'GMBSV_GMSV_FK' },
        build             => { is => 'Genome::Model::Build', id_by => 'build_id', constraint_name => 'GMBSV_GMB_FK' },
        allele_frequency  => { is => 'NUMBER', len => 12 },
        breakdancer_score => { is => 'NUMBER', len => 12 },
        num_reads         => { is => 'NUMBER', len => 12 },
        num_reads_lib     => { is => 'VARCHAR2', len => 255 },
        run_param         => { is => 'VARCHAR2', len => 255, is_optional => 1 },
        somatic_status    => { is => 'VARCHAR2', len => 20 },
        version           => { is => 'VARCHAR2', len => 50 },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;
