package Genome::Model::SomaticValidation::Command::ManualResult;

use strict;
use warnings;

use Genome;
use Cwd;
use Genome::Model::Tools::DetectVariants2::Utilities qw(
    final_result_for_variant_type
);

class Genome::Model::SomaticValidation::Command::ManualResult {
    is => 'Command::V2',
    has_input => [
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
        source_build => {
            id_by => 'source_build_id',
            is => 'Genome::Model::Build',
            doc => 'The build containing the data on which these variant calls were based',
        },
        reference_sequence_build => {
            id_by => 'reference_sequence_build_id',
            is => 'Genome::Model::Build',
            doc => 'The reference build used for these variants - only required if source_build is not given',
        },
        format => {
            is => 'Text',
            doc => 'The format of the variant list (e.g. "bed", "samtools", "breakdancer")',
            default_value => 'bed',
        },
        analysis_project => {
            is => 'Genome::Config::AnalysisProject',
            doc => 'The Analysis Project for which the manual result is needed',
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

    $self->variant_file(Cwd::abs_path($self->variant_file));

    if (!defined($self->source_build) && !defined($self->reference_sequence_build)){
        die("Must supply either a source build or a reference build");
    }


    my %params = (
        variant_type => $self->variant_type,
        original_file_path => $self->variant_file,
        description => $self->description,
        test_name => Genome::Config::get('software_result_test_name') || undef,
        format => $self->format,
    );

    my $manual_result;

    if(defined($self->source_build)){
        my $source_build = $self->source_build;
        my $previous_result = final_result_for_variant_type([$source_build->results], $self->variant_type . 's');
        $params{previous_result_id} = ($previous_result? $previous_result->id : undef);
        $params{source_build_id} = $source_build->id;

        # allow tumor only or normal only models.
        if ($source_build->model->experimental_subject) {
            $params{sample_id} = $source_build->model->experimental_subject->id;
        }
        if ($source_build->model->control_subject) {
            $params{control_sample_id} = $source_build->model->control_subject->id;
        }
        $params{reference_build_id} = $source_build->reference_sequence_build->id;

        $params{users} = {
            requestor => $source_build,
            sponsor   => ($self->analysis_project // Genome::Sys->current_user),
        };
        $manual_result = Genome::Model::Tools::DetectVariants2::Result::Manual->get_or_create(%params);
    } else {
        $params{reference_build_id} = $self->reference_sequence_build->id;
        $manual_result = Genome::Model::Tools::DetectVariants2::Result::Manual->create(%params);
    }


    unless($manual_result) {
        die $self->error_message('Failed to generate new result for data.');
    }

    $self->manual_result($manual_result);
    $self->status_message('Created a manual result.  The ID for this new result is: ' . $manual_result->id);

    return 1;
}

1;
