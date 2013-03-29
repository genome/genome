package Genome::Model::VariantSource;

use strict;
use warnings;

use Genome;
class Genome::Model::VariantSource {
    table_name => 'GENOME_MODEL_VARIANT_SOURCE',
    er_role => 'bridge',
    type_name => 'genome model variant source',
    id_by => [
        build_id => {
            is => 'NUMBER',
            len => 10,
        },
        source_id => {
            is => 'NUMBER',
            len => 10,
        },
        var_id => {
            is => 'NUMBER',
            len => 10,
        },
    ],
    has => [
        genome_model_build => {
            is => 'Genome::Model::Build',
            id_by => 'build_id',
            constraint_name => 'FK_GMB_BI',
        },
        genome_model_variant_var => {
            is => 'Genome::Model::Variant',
            id_by => 'var_id',
            constraint_name => 'FK_GMV_VI',
        },
        genome_variant_source => {
            is => 'Genome::VariantSource',
            constraint_name => 'FK_GVS_SI',
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;
