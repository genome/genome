# Review: gsanders - Some sort of enforcement needs to be made on validation_type and validation_result
# The current practice is to use "Solexa" and "Official" for validation_type and "S" "WT" "P" etc for validation_result.
# Also we have run into problems recently with "S " (with a space or a newline) being uploaded to db. Implement valid_values probably.

package Genome::Model::VariantValidation;

use strict;
use warnings;

use Genome;
class Genome::Model::VariantValidation {
    type_name => 'genome model variantvalidation',
    table_name => 'GENOME_MODEL_VARIANTVALIDATION',
    id_by => [
        model_id        => { is => 'NUMBER', len => 11, implied_by => 'model' },
        validation_type => { is => 'VARCHAR2', 
                             len => 255, 
                             valid_values => ["Illumina", "3730", "Capture", "454", "Official"],
                           },
        variant_id      => { is => 'NUMBER', len => 10, implied_by => 'variant' },
    ],
    has => [
        model             => { is => 'Genome::Model', id_by => 'model_id', constraint_name => 'GMV_GM_FK' },
        variant           => { is => 'Genome::Model::Variant', id_by => 'variant_id', constraint_name => 'GMV_GMV_FK' },
        validation_result => { is => 'VARCHAR2', len => 10,
                                valid_values => ["P", "WT", "G", "S", "LOH", "NC"],
                             },
        comments          => { is => 'VARCHAR2', len => 255, is_optional => 1 },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;
