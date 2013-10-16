package Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePath;

use strict;
use warnings;

use Genome;

require File::Basename;

class Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePath { 
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
        destination_path => {
            calculate_from => [qw/ working_directory source_base_name /],
            calculate => q( return $working_directory.'/'.$source_base_name; ),
            doc => 'Final destination path.',
        }, 
    ],
    has_optional_calculated => [
        source_base_name => {
            calculate_from => [qw/ source_path /],
            calculate => q( return File::Basename::basename($source_path); ),
        },
    ],
    has_constant_calculated => [
        helpers => {
            calculate => q( Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get; ),
        },
    ],
};

sub execute {
    my $self = shift;

    my $retrieve_source_path = $self->_retrieve_source_path;
    return if not $retrieve_source_path;

    my $retrieve_source_md5 = $self->_retrieve_source_md5;
    return if not $retrieve_source_md5;

    return 1;
}

sub _retrieve_path {
    my ($self, $source_path, $destination_path) = @_;

    $self->status_message('Source: '.$source_path);
    $self->status_message('Destination: '.$destination_path);
    my $copy_ok = $self->helpers->copy_file($source_path, $destination_path);
    return $copy_ok
}

sub _retrieve_source_path {
    my $self = shift;
    $self->status_message('Retrieve source path...');

    my $retrieve_ok = $self->_retrieve_path($self->source_path, $self->destination_path);
    return if not $retrieve_ok;
    
    $self->status_message('Retrieve source path...done');
    return 1;
}

sub _retrieve_source_md5 {
    my $self = shift;
    $self->status_message('Retrieve source MD5 path...');

    my $md5_path = $self->source_path.'.md5';
    my $md5_size = $self->helpers->file_size($md5_path);
    if ( not $md5_size ) {
        $self->status_message('Source MD5 path does not exist...skip');
        return 1;
    }

    my $retrieve_ok = $self->_retrieve_path($md5_path, $self->destination_path.'.md5-orig');
    return if not $retrieve_ok;

    $self->status_message('Retrieve source MD5 path...done');
    return 1;
}

1;

