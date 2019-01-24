package Genome::FeatureList::Command::DumpIntervalList;

use strict;
use warnings;

use Genome;

class Genome::FeatureList::Command::DumpIntervalList {
    is => 'Command::V2',
    has_input => {
        feature_list => {
            is => 'Genome::FeatureList',
            doc => 'the feature list to convert to interval list',
        },
        reference_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            doc => 'the reference for which to produce the interval list',
        },
        track_name => {
            is => 'Text',
            doc => 'For multi-tracked BED files, which track to use',
            valid_values => Genome::FeatureList::IntervalList->__meta__->property(property_name => 'track_name')->valid_values,
        },
        merge => {
            is => 'Boolean',
            doc => 'whether to merge adjacent regions of the BED file',
            default => 1,
        },
        preserve_region_names => {
            is => 'Boolean',
            doc => 'whether to retain the region names from the feature-list (or replace them with "short" names)',
            default => 0,
        },
    },
    has_optional_output => {
        output_path => {
            is => 'Path',
            doc => 'where to write the interval_list (default: STDOUT)',
        },
    },
};

sub execute {
    my $self = shift;

    my $result_users = {
        sponsor => Genome::Config::AnalysisProject->system_analysis_project,
        requestor => $self->reference_build,
    };

    my $sr = Genome::FeatureList::IntervalList->get_or_create(
        feature_list => $self->feature_list,
        reference_build => $self->reference_build,
        track_name => $self->track_name,
        merge => $self->merge,
        preserve_region_names => $self->preserve_region_names,
        users => $result_users,
    );

    unless($sr) {
        $self->fatal_message('Failed to get or create interval list.');
    }
    my $file = $sr->interval_list;

    if ($self->output_path) {
        Genome::Sys->copy_file($file, $self->output_path);
    } else {
        Genome::Sys->shellcmd(
            cmd => ['cat', $file],
            input_files => [$file],
        );
    }

    return 1;
}

1;
