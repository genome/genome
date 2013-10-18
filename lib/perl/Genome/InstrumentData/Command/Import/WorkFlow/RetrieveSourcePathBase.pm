package Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathBase;

use strict;
use warnings;

use Genome;

require File::Basename;

class Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathBase { 
    is => 'Command::V2',
    has_input => [
        working_directory => {
            is => 'Text',
            doc => 'Detination directory for source path.',
        },
    ],
    has_constant_calculated => [
        helpers => {
            calculate => q( Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get; ),
        },
    ],
};

sub destination_path_for_source_path {
    my ($self, $source_path) = @_;
    my $source_base_name = File::Basename::basename($source_path);
    return $self->working_directory.'/'.$source_base_name;
}

sub _retrieve_source_path {
    my ($self, $source_path) = @_;
    $self->status_message('Retrieve source path...');

    Carp::confess('No source path to retrieve!') if not $source_path;

    my $retrieve_ok = $self->helpers->retrieve_path($source_path, $self->destination_path_for_source_path($source_path));
    return if not $retrieve_ok;
    
    $self->status_message('Retrieve source path...done');
    return 1;
}

sub _retrieve_source_md5 {
    my ($self, $source_path) = @_;
    $self->status_message('Retrieve source MD5 path...');

    Carp::confess('No source path to retrieve MD5!') if not $source_path;

    my $md5_path = $source_path.'.md5';
    my $md5_size = $self->helpers->file_size($md5_path);
    if ( not $md5_size ) {
        $self->status_message('Source MD5 path does not exist...skip');
        return 1;
    }

    my $retrieve_ok = $self->helpers->retrieve_path($md5_path, $self->destination_path_for_source_path($source_path).'.md5-orig');
    return if not $retrieve_ok;

    $self->status_message('Retrieve source MD5 path...done');
    return 1;
}

1;

