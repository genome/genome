package Genome::InstrumentData::Command::Import::WorkFlow::VerifyNotImported;

use strict;
use warnings;

use Genome;

require File::Basename;
require Genome::InstrumentData::Command::Import::Inputs::SourceFile;

class Genome::InstrumentData::Command::Import::WorkFlow::VerifyNotImported { 
    is => [qw/ Command::V2 Genome::Model::Tools::Picard::WithDownsampleRatio /],
    has_input => {
        working_directory => {
            is => 'Text',
            doc => 'Detination directory for source path.',
        },
        source_path => {
            is => 'Text',
            is_output => 1,
            doc => 'Source path of sequences to get.',
            is_output => 1,
        },
    },
    has_output => {
        source_md5 => {
            is => 'Text',
            doc => 'Source MD5.',
        }, 
    },
    has_constant_calculated => {
        helpers => {
            calculate => q( Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get; ),
        },
        source_file => {
            is_constant => 1,
            calculate_from => [qw/ source_path /],
            calculate => q| return Genome::InstrumentData::Command::Import::Inputs::SourceFile->create(path => $source_path); |,
        },
        output_file => {
            is_constant => 1,
            calculate_from => [qw/ working_directory source_path /],
            calculate => q| 
                return Genome::InstrumentData::Command::Import::Inputs::SourceFile->create(
                    path => File::Spec->join($working_directory, File::Basename::basename($source_path)),
                ); 
            |,
        },
    },
};

sub execute {
    my $self = shift;

    # Copy MD5 path for source file to working directory, if exists
    my $md5_path = $self->output_file->md5_path;
    if ( $self->source_file->md5_path_size ) {
        Genome::Sys->copy_file($self->source_file->md5_path, $md5_path);
    }
    # Otherwise run it, saving to working directory
    else {
        $self->helpers->run_md5($self->output_file->path, $md5_path);
    }
    return if not -s $md5_path;

    # Load MD5
    my $md5 = $self->helpers->load_md5($md5_path);

    # Verify md5 not previouly imported
    my %downsample_ratio_param;
    $downsample_ratio_param{downsample_ratio} = $self->downsample_ratio if defined $self->downsample_ratio;
    my $previously_imported = $self->helpers->were_original_path_md5s_previously_imported(
        md5s => [ $md5 ],
        %downsample_ratio_param,
    );
    return if $previously_imported;

    $self->source_md5($md5);
    return 1;
}

1;

