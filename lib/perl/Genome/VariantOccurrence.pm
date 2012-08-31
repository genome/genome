# review gsanders jlolofie
# Talked to charris... it seems there are ~100 rows in this table. This is old validation data. What needs to happen is that this data
# should be imported into the new Genome::Model::Variant, Genome::Model::BuildVariant, Genome::Model::VariantValidation schema.
# This will take a little manual work as we do not have the full annotation for these variants handy so we need to re-annotate these lines 
# and format them to be uploaded into the new table set. This will take some manual munging but nothing too bad
# Once this is done we can drop this table and this class.

package Genome::VariantOccurrence;

use strict;
use warnings;

use Genome;
class Genome::VariantOccurrence {
    type_name => 'variant occurrence',
    table_name => 'VARIANT_OCCURRENCE',
    id_by => [
        build_id   => { is => 'NUMBER', len => 10, implied_by => 'genome_model_build' },
        chromosome => { is => 'VARCHAR2', len => 50 },
        start_pos  => { is => 'NUMBER', len => 10 },
        stop_pos   => { is => 'NUMBER', len => 10 },
    ],
    has => [
        gene               => { is => 'VARCHAR2', len => 50, is_optional => 1 },
        genome_model_build => { is => 'Genome::Model::Build', id_by => 'build_id', constraint_name => 'VO_GMB_FK' },
        reference_allele   => { is => 'VARCHAR2', len => 50 },
        somatic_status     => { is => 'VARCHAR2', len => 1 },
        variant_allele     => { is => 'VARCHAR2', len => 50 },
        type               => { is => 'VARCHAR2', len => 3, is_optional => 1 },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;
