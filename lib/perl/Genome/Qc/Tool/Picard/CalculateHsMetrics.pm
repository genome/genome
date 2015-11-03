package Genome::Qc::Tool::Picard::CalculateHsMetrics;

use strict;
use warnings;
use Genome;

use List::MoreUtils qw(uniq);

class Genome::Qc::Tool::Picard::CalculateHsMetrics {
    is => 'Genome::Qc::Tool::Picard',
};

sub supports_streaming {
    return 1;
}

sub qc_metrics_file_accessor {
    return 'output_file';
}

sub gmt_class {
    return 'Genome::Model::Tools::Picard::CalculateHsMetrics';
}

sub bait_intervals {
    my $self = shift;

    my $translated_list = $self->translate_feature_list('tiled_region');
    return $self->create_interval_file($translated_list);
}

sub target_intervals {
    my $self = shift;

    my $translated_list = $self->translate_feature_list('target_region');
    return $self->create_interval_file($translated_list);
}

sub translate_feature_list {
    my ($self, $track_name) = @_;

    my $feature_list = $self->_feature_list_for_alignment_result($self->alignment_result);

    my $translated_list = Genome::Sys->create_temp_file_path;
    my %params = (
        feature_list => $feature_list,
        track_name => $track_name,
        output_path => $translated_list,
        merge => 0,
    );
    unless ($self->reference_build eq $feature_list->reference) {
        $params{alternate_reference} = $self->reference_build;
    }
    my $translated_list_cmd = Genome::FeatureList::Command::DumpMergedList->create(%params);
    $translated_list_cmd->execute;

    return $translated_list;
}

sub _feature_list_for_alignment_result {
    my $self = shift;
    my $alignment_result = shift;

    my @instrument_data = $alignment_result->instrument_data;
    my @trsn = uniq map { $_->target_region_set_name } @instrument_data;
    unless(@trsn) {
        $self->fatal_message(
            'Alignment result %s does not have an associated target region set.',
            $alignment_result->__display_name__
        );
    } elsif(@trsn > 1) {
        $self->fatal_message(
            'Alignment result %s contains instrument data with different target region sets.',
            $alignment_result->__display_name__
        );
    }

    return $instrument_data[0]->target_region_set;
}

sub create_interval_file {
    my ($self, $translated_list) = @_;

    my $interval_file = Genome::Sys->create_temp_file_path;
    my $interval_file_cmd = Genome::Model::Tools::Bed::ToIntervals->create(
        bed_file => $translated_list,
        interval_file => $interval_file,
        seqdict_file => $self->reference_build->get_sequence_dictionary('sam'),
    );
    $interval_file_cmd->execute;

    return $interval_file;
}

1;
