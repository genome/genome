package Genome::InstrumentData::Command::Import::WorkFlow::VerifyNotImported;

use strict;
use warnings;

use Genome;

require File::Basename;

class Genome::InstrumentData::Command::Import::WorkFlow::VerifyNotImported { 
    is => [qw/ Command::V2 Genome::Model::Tools::Picard::WithDownsampleRatio /],
    has_input => {
        working_directory => {
            is => 'Text',
            doc => 'Detination directory for source path.',
        },
        source_path => {
            is => 'Text',
            doc => 'Source path of sequences to get.',
        },
    },
    has_output => {
        source_md5 => {
            is => 'Text',
            doc => 'Source MD5.',
        }, 
    },
    has_optional_calculated => {
        source_md5_path => {
            calculate_from => [qw/ source_path /],
            calculate => q{ return $self->helpers->md5_path_for($source_path); },
            doc => 'Source MD5 path.',
        }, 
        original_md5_path => {
            calculate_from => [qw/ source_path /],
            calculate => q{ return $self->helpers->original_md5_path_for($source_path); },
            doc => 'Original source MD5 path.',
        }, 
    },
    has_transient_optional => {
        original_md5 => { is => 'Text', },
    },
    has_constant_calculated => {
        helpers => {
            calculate => q( Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get; ),
        },
    },
};

sub execute {
    my $self = shift;

    # Load original md5, if exists
    my $load_original_md5 = $self->_load_original_md5;
    return if not $load_original_md5;

    my %downsample_ratio_param;
    $downsample_ratio_param{downsample_ratio} = $self->downsample_ratio if defined $self->downsample_ratio;

    # Verify original md5 [if exists] not previously imported
    my $original_md5 = $self->original_md5;
    if ( $original_md5 ) { # check if previously imported
        my $previously_imported = $self->helpers->were_original_path_md5s_previously_imported(
            md5s => [ $original_md5 ],
            %downsample_ratio_param,
        );
        return if $previously_imported;
    }

    # Run md5 on source file
    my $md5 = $self->helpers->load_or_run_md5($self->source_path, $self->source_md5_path);
    return if not $md5;

    # Verify original md5 [if exists] matches
    if ( $original_md5 and $md5 ne $original_md5 ) {
        $self->error_message("Original and generated MD5s do not match! $original_md5 vs. $md5");
        return;
    }

    # Verify md5 not previouly imported
    my $previously_imported = $self->helpers->were_original_path_md5s_previously_imported(
        md5s => [ $md5 ],
        %downsample_ratio_param,
    );
    return if $previously_imported;

    $self->source_md5($md5);
    return 1;
}

sub _load_original_md5 {
    my $self = shift;
    $self->debug_message('Load original MD5...');

    my $original_md5_path = $self->helpers->original_md5_path_for($self->source_path);
    my $original_md5_path_size = $self->helpers->file_size($original_md5_path);
    if ( not $original_md5_path_size ) {
        $self->debug_message('No original MD5...skip');
        return 1;
    }

    my $original_md5 = $self->helpers->load_md5($original_md5_path);
    if ( not $original_md5 ) {
        $self->error_message('Failed to load original MD5!');
        return;
    }
    $self->original_md5($original_md5);
    $self->debug_message("Original MD5: $original_md5");

    $self->debug_message('Load original MD5...done');
    return 1;
}

1;

