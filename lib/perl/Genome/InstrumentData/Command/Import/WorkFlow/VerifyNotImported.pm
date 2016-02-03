package Genome::InstrumentData::Command::Import::WorkFlow::VerifyNotImported;

use strict;
use warnings;

require File::Basename;

use Genome;
require Genome::InstrumentData::Command::Import::Inputs::SourceFile;

class Genome::InstrumentData::Command::Import::WorkFlow::VerifyNotImported { 
    is => [qw/ Command::V2 Genome::Model::Tools::Picard::WithDownsampleRatio /],
    roles => [qw/ Genome::InstrumentData::Command::Import::WorkFlow::Role::WithWorkingDirectory /],
    has_input => {
        source_paths => {
            is => 'Text',
            is_many => 1,
            doc => 'Source path of sequences to get.',
        },
    },
    has_output => {
        output_path => {
            is => 'Text',
            calculate_from => [qw/ _output_paths /],
            calculate => q| return $_output_paths->[0]; |,
        },
        output_paths => {
            is => 'Text',
            is_many => 1,
            calculate_from => [qw/ _output_paths /],
            calculate => q| return @$_output_paths; |,
        },
        source_md5s => {
            is => 'Text',
            is_many => 1,
            doc => 'Source MD5.',
        }, 
    },
    has_constant_calculated => {
        helpers => {
            calculate => q( Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get; ),
        },
    },
    has_optional_transient => {
        _output_paths => { is => 'ARRAY', },
    },
};

sub source_file_for_path {
    my ($self, $path) = @_;
    return Genome::InstrumentData::Command::Import::Inputs::SourceFile->create(path => $path);
}

sub output_file_for_path {
    my ($self, $path) = @_;
    return Genome::InstrumentData::Command::Import::Inputs::SourceFile->create(
        path => File::Spec->join($self->working_directory, File::Basename::basename($path)),
    ); 
}

sub execute {
    my $self = shift;

    my (@output_paths, @md5s);
    for my $source_path ( $self->source_paths ) {
        my $source_file = $self->source_file_for_path($source_path);
        my $output_file = $self->output_file_for_path($source_path);

        # Link source path to working_directory
        $self->fatal_message('Output file path exists! Cannot symlink to it! Are the source files in the working directory?') if -e $output_file->path;
        Genome::Sys->create_symlink($source_path, $output_file->path) if not -e $output_file->path;
        push @output_paths, $output_file->path;

        # Copy MD5 path for source file to working directory, if exists
        my $md5_path = $output_file->md5_path;
        if ( $source_file->md5_path_size ) {
            Genome::Sys->copy_file($source_file->md5_path, $md5_path);
        }
        # Otherwise run it, saving to working directory
        else {
            $self->helpers->run_md5($output_file->path, $md5_path);
        }
        return if not -s $md5_path;

        # Load MD5
        my $md5 = $self->helpers->load_md5($md5_path);
        push @md5s, $md5;

        # Verify md5 not previouly imported
        my %downsample_ratio_param;
        $downsample_ratio_param{downsample_ratio} = $self->downsample_ratio if defined $self->downsample_ratio;
        my $previously_imported = $self->helpers->were_original_path_md5s_previously_imported(
            md5s => [ $md5 ],
            %downsample_ratio_param,
        );
        return if $previously_imported;
    }

    $self->_output_paths(\@output_paths);
    $self->source_md5s(\@md5s);
    return 1;
}

1;

