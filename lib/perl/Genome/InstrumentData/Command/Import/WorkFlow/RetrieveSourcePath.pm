package Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePath;

use strict;
use warnings;

use Genome;

require File::Basename;
require File::Spec;

class Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePath { 
    is => 'Command::V2',
    is_abstract => 1,
    has_input => {
        source_path => {
            is => 'Text',
            doc => 'Source path of sequences to get.',
        },
       working_directory => {
            is => 'Text',
            doc => 'Detination directory for source path.',
        },
    },
    has_output => {
        destination_path => {
            calculate_from => [qw/ working_directory source_path /],
            calculate => q| return $self->working_directory.'/'.File::Basename::basename($source_path); |,
            doc => 'Final destination path.',
        }, 
        destination_md5_path => {
            calculate_from => [qw/ destination_path /],
            calculate => q| return Genome::InstrumentData::Command::Import::WorkFlow::Helpers->original_md5_path_for($destination_path); |,
            doc => 'Final destination path.',
        }, 
    },
    has_constant_calculated => {
        helpers => {
            calculate => q( Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get; ),
        },
    },
};

sub execute {
    my $self = shift;

    my $retrieve_source_path = $self->_retrieve_source_path;
    return if not $retrieve_source_path;

    my $retrieve_source_md5 = $self->_retrieve_source_md5_path;
    return if not $retrieve_source_md5;

    return 1;
}

sub _retrieve_source_path {
    my $self = shift;
    return $self->_retrieve_path($self->source_path, $self->destination_path);
}

sub _retrieve_source_md5_path {
    my $self = shift;

    my $md5_path = $self->helpers->md5_path_for($self->source_path);
    my $md5_size = $self->helpers->file_size($md5_path);
    if ( not $md5_size ) {
        $self->debug_message('Source MD5 path does not exist...skip');
        return 1;
    }

    $self->debug_message('Retrieve source MD5 path...');

    my $original_md5_path = $self->helpers->original_md5_path_for($self->destination_path);
    my $retrieve_ok = $self->_retrieve_path($md5_path, $original_md5_path);
    return if not $retrieve_ok;

    $self->debug_message('Retrieve source MD5 path...done');
    return 1;
}

1;

