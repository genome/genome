package Genome::FeatureList::IntervalList;

use strict;
use warnings;

use Genome;

use File::Spec;

class Genome::FeatureList::IntervalList {
    is => ['Genome::SoftwareResult::StageableSimple', 'Genome::SoftwareResult::WithNestedResults'],
    has_input => {
        feature_list => {
            is => 'Genome::FeatureList',
            doc => 'the feature list to convert to interval list',
        },
        reference_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            doc => 'the reference in which to produce the interval list',
        },
    },
    has_param => [
        track_name => {
            is => 'Text',
            doc => 'For multi-tracked BED files, which track to use',
            valid_values => Genome::FeatureList::Command::DumpMergedList->__meta__->property(property_name => 'track_name')->valid_values,
        },
        merge => {
            is => 'Boolean',
            doc => 'whether to merge adjacent regions of the BED file',
            default => 1,
        },
        short_name => {
            is => 'Boolean',
            doc => 'whether to replace region names with short names',
            default => 1,
        },
    ],
};

sub interval_list {
    my $self = shift;

    return File::Spec->join($self->output_dir, $self->_file_name);
}

sub _file_name {
    my $self = shift;
    my $fl = $self->feature_list;

    return $fl->id . '.interval_list';
}

sub _run {
    my $self = shift;

    my $dumped_bed_path = Genome::Sys->create_temp_file_path;

    my @alt_reference;
    if ($self->reference_build ne $self->feature_list->reference) {
        @alt_reference = (alternate_reference => $self->reference_build);
    }

    my $bed_dumper = Genome::FeatureList::Command::DumpMergedList->create(
        output_path => $dumped_bed_path,
        feature_list => $self->feature_list,
        track_name => $self->track_name,
        merge => $self->merge,
        short_name => $self->short_name,
        result_users => $self->_user_data_for_nested_results,
        @alt_reference,
    );
    unless ($bed_dumper and $bed_dumper->execute) {
        $self->fatal_message('Failed to get dumped BED file for conversion');
    }

    my $interval_list_path = File::Spec->join($self->temp_staging_directory, $self->_file_name);

    my $converter = Genome::Model::Tools::Picard::BedToIntervalList->create(
        input => $dumped_bed_path,
        output => $interval_list_path,
        sequence_dictionary => $self->reference_build->sequence_dictionary_path('sam'),
        use_version => '1.124',
    );
    unless($converter and $converter->execute) {
        $self->fatal_message('Failed to convert BED to interval list');
    }

    return 1;
}

sub resolve_allocation_subdirectory {
    my $self = shift;

    return File::Spec->join('model_data', 'interval-list', $self->id);
}

sub resolve_allocation_disk_group_name {
    Genome::Config::get('disk_group_references');
}

1;
