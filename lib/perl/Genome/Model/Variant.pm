# Review: gsanders - Validation status needs to be removed... this is all tracked (including the "official" call) in variantvalidation

package Genome::Model::Variant;

use strict;
use warnings;

use Genome;
class Genome::Model::Variant {
    type_name => 'genome model variant',
    table_name => 'GENOME_MODEL_VARIANT',
    id_by => [
        variant_id => { is => 'NUMBER', len => 10 },
    ],
    has => [
        chromosome         => { is => 'VARCHAR2', len => 255 },
        start_pos          => { is => 'NUMBER', len => 12 },
        stop_pos           => { is => 'NUMBER', len => 12 },
        reference_allele   => { is => 'VARCHAR2', len => 255, is_optional => 1 },
        variant_allele     => { is => 'VARCHAR2', len => 255, is_optional => 1 },
        amino_acid_change  => { is => 'VARCHAR2', len => 255, is_optional => 1 },
        c_position         => { is => 'VARCHAR2', len => 255, is_optional => 1 },
        domain             => { is => 'VARCHAR2', len => 255, is_optional => 1 },
        gene_name          => { is => 'VARCHAR2', len => 255, is_optional => 1 },
        strand             => { is => 'VARCHAR2', len => 255, is_optional => 1 },
        transcript_name    => { is => 'VARCHAR2', len => 255, is_optional => 1 },
        transcript_source  => { is => 'VARCHAR2', len => 255, is_optional => 1 },
        transcript_status  => { is => 'VARCHAR2', len => 255, is_optional => 1 },
        transcript_version => { is => 'VARCHAR2', len => 255, is_optional => 1 },
        trv_type           => { is => 'VARCHAR2', len => 255, is_optional => 1 },
        ucsc_cons          => { is => 'VARCHAR2', len => 255, is_optional => 1 },
        validation_links   => { is => 'Genome::Model::VariantValidation', reverse_as => 'variant', is_many => 1 },
    ],
    unique_constraints => [
        { properties => [qw/chromosome reference_allele start_pos stop_pos variant_allele/], sql => 'GMV_UK' },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

# This method uses the build id passed in to get the "official" validation for this variant on the model id to which the build belongs
sub get_official_validation_for_build {
    my $self = shift;
    my $build_id = shift;

    my $build = Genome::Model::Build->get($build_id);
    unless($build) {
        $self->error_message("Could not get a build for build id: " . $self->build_id . ". Please use a valid build id.");
        return;
    }
    my $model = $build->model;
    unless($model) {
        $self->error_message("Could not get a model for model id: " . $build->model_id . ". Please use a build with a valid model id.");
        return;
    }

    my @validations = $self->validation_links;
    my $official_validation;
    for my $validation (@validations) {
        if (($validation->model_id == $model->genome_model_id)&&($validation->validation_type eq 'Official')) {
            # There should only be one official validation
            if ($official_validation) {
                $self->error_message("More than one 'Official' validation found for model_id " . $build->model_id);
                return;
            } else {
                $official_validation = $validation;
            }
        }
    }

    unless ($official_validation) {
        $self->warning_message("No 'Official' validation found for model_id " . $build->model_id);
        return;
    }

    return $official_validation;
}

1;
