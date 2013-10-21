package Genome::InstrumentData::Command::Import::WorkFlow::VerifyMd5;

use strict;
use warnings;

use Genome;

require File::Basename;

class Genome::InstrumentData::Command::Import::WorkFlow::VerifyMd5 { 
    is => 'Command::V2',
    has_input => [
        working_directory => {
            is => 'Text',
            doc => 'Detination directory for source path.',
        },
        source_path => {
            is => 'Text',
            doc => 'Source path of sequences to get.',
        },
    ],
    has_output => [
        source_md5_path => {
            calculate_from => [qw/ working_directory source_path_base_name /],
            calculate => q( return $working_directory.'/'.$source_path_base_name.'.md5'; ),
            doc => 'Source MD5 path.',
        }, 
    ],
    has_optional_calculated => [
        source_path_base_name => {
            calculate_from => [qw/ source_path /],
            calculate => q( return File::Basename::basename($source_path); ),
        },
        original_md5_path => {
            calculate_from => [qw/ working_directory source_path_base_name /],
            calculate => q( return $working_directory.'/'.$source_path_base_name.'.orig-md5'; ),
        },
    ],
    has_transient_optional => [
        original_md5 => { is => 'Text', },
    ],
    has_constant_calculated => [
        helpers => {
            calculate => q( Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get; ),
        },
    ],
};

sub execute {
    my $self = shift;

    my $load_original_md5 = $self->_load_original_md5;
    return if not $load_original_md5;

    my $original_md5 = $self->original_md5;
    my $md5;
    if ( $original_md5 ) {
        my @instrument_data_attr = Genome::InstrumentDataAttribute->get(
            attribute_label => 'original_data_path_md5',
            attribute_value => $original_md5,
        );
        if ( @instrument_data_attr ) {
            $self->error_message(
                "Instrument data was previously imported! Found existing instrument data with MD5 ($original_md5): ".
                join(' ', map { $_->instrument_data_id } @instrument_data_attr)
            );
            return;
        }

        $md5 = $self->helpers->load_or_run_md5($self->source_path, $self->source_md5_path);
        return if not $md5;

        if ( $md5 ne $original_md5 ) {
            $self->error_message("Original and generated MD5s do not match! $original_md5 vs. $md5");
            return;
        }
    }
    else {
        $md5 = $self->helpers->load_or_run_md5($self->source_path, $self->source_md5_path);
        return if not $md5;
    }

    return 1;
}

sub _load_original_md5 {
    my $self = shift;
    $self->status_message('Load original MD5...');

    my $original_md5_path = $self->original_md5_path;
    my $original_md5_path_size = $self->helpers->file_size($original_md5_path);
    if ( not $original_md5_path_size ) {
        $self->status_message('No original MD5...skip');
        return 1;
    }

    my $original_md5 = $self->helpers->load_md5($original_md5_path);
    if ( not $original_md5 ) {
        $self->error_message('Failed to load original MD5!');
        return;
    }
    $self->original_md5($original_md5);
    $self->status_message("Original MD5: $original_md5");

    $self->status_message('Load original MD5...done');
    return 1;
}

1;

