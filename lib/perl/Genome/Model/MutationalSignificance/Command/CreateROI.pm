package Genome::Model::MutationalSignificance::Command::CreateROI;

use strict;
use warnings;

use Genome;

class Genome::Model::MutationalSignificance::Command::CreateROI {
    is => ['Command::V2'],
    has_input => [
        annotation_build => {
            is => 'Genome::Model::Build::ImportedAnnotation'
        },
        excluded_reference_sequence_patterns => {
            is => 'Text',
            is_many => 1,
            is_optional => 1,
            doc => "Exclude transcripts on these reference sequences",
        },
        included_feature_type_patterns => {
            is => 'Text',
            is_many => 1,
            is_optional => 1,
            doc => 'Include only entries that match one of these patterns',
        },
        condense_feature_name => {
            is => 'Boolean',
            doc => 'Use only gene name as feature name',
            default_value => 1,
        },
        flank_size => {
            is => 'Integer',
            doc => 'Add this number of base pairs on each side of the feature', #to do: check this
            default_value => 0,
        },
        one_based => {
            is => 'Boolean',
            default_value => 1,
        },
    ],
    has_output => [
        roi_path => {
            is => 'String',
        },
    ],
};

sub execute {
    my $self = shift;

    my @params;

    if ($self->excluded_reference_sequence_patterns) {
        push @params, excluded_reference_sequence_patterns => [$self->excluded_reference_sequence_patterns];
    }
    if ($self->included_feature_type_patterns) {
        push @params, included_feature_type_patterns => [$self->included_feature_type_patterns];
    }
    if ($self->condense_feature_name) {
        push @params, condense_feature_name => $self->condense_feature_name;
    }
    if ($self->flank_size && $self->flank_size > 0) {
        push @params, flank_size => $self->flank_size;
    }

    my $feature_list = $self->annotation_build->get_or_create_roi_bed(@params);

    unless ($feature_list) {
        $self->error_message('ROI file not available from annotation build '.$self->annotation_build->id);
        return;
    }

    $self->roi_path($feature_list->file_path);

    $self->status_message('Using ROI file: '.$self->roi_path);
    return 1;
}

1;
