package Genome::Model::SomaticValidation::Command::ManualResult;

use strict;
use warnings;

use Genome;
use Cwd;

class Genome::Model::SomaticValidation::Command::ManualResult {
    is => 'Command::V2',
    has_input => [
        source_build => {
            id_by => 'source_build_id',
            is => 'Genome::Model::Build',
            doc => 'The build on which these variants are based',
        },
        variant_file => {
            is => 'FilePath',
            doc => 'Path to the file of variants',
        },
        variant_type => {
            is => 'Text',
            doc => 'The type of variants in this result',
            valid_values => ['snv', 'indel', 'sv', 'cnv'],
        },
        description => {
            is => 'Text',
            doc => 'General description of the list',
        },
    ],
    has_optional_input => [
        format => {
            is => 'Text',
            doc => 'The format of the variant list (e.g. "bed", "samtools", "breakdancer")',
            default_value => 'bed',
        },
    ],
    has_transient_optional_output => [
        manual_result => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Manual',
            doc => 'The SoftwareResult record created for the supplied data',
        },
    ],
    doc => 'log a reduced list of variants selected for validation',

};

sub sub_command_category { 'analyst tools' }

sub execute {
    my $self = shift;

    my $source_build = $self->source_build;
    my $previous_result;
    if($source_build->can('final_result_for_variant_type')){
        $previous_result = $source_build->final_result_for_variant_type($self->variant_type .'s');
    }

    my $tumor_model = $source_build->can('tumor_model')? $source_build->tumor_model : $source_build->model->tumor_model;
    my $normal_model = $source_build->can('normal_model')? $source_build->normal_model : $source_build->model->normal_model;

    $self->variant_file(Cwd::abs_path($self->variant_file));

    my $manual_result = Genome::Model::Tools::DetectVariants2::Result::Manual->get_or_create(
        variant_type => $self->variant_type,
        sample_id => $tumor_model->subject->id,
        control_sample_id => $normal_model->subject->id,
        reference_build_id => $source_build->reference_sequence_build->id,
        original_file_path => $self->variant_file,
        description => $self->description,
        format => $self->format,
        previous_result_id => ($previous_result? $previous_result->id : undef),
        test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef,
        source_build_id => $source_build->id,
    );

    unless($manual_result) {
        die $self->error_message('Failed to generate new result for data.');
    }

    $self->manual_result($manual_result);
    $self->status_message('Created a manual result.  The ID for this new result is: ' . $manual_result->id);

    return 1;
}

1;
