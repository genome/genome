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
};

sub execute {
    my $self = shift;
    $self->status_message('Retrieve source path...');

    # TODO
    # add more methods to get
    # md5 verification

    my $retrieve_source_path = $self->_retrieve_source_path;
    return if not $retrieve_source_path;

    $self->status_message('Retrieve source path...done');
    return 1;
}

sub _retrieve_source_path {
    my ($self, $s) = @_;

    my $source_path = $self->source_path;
    $self->status_message('Source path: '.$source_path);
    my $destination_path = $self->destination_path;
    $self->status_message('Destination path: '.$destination_path);

    my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;
    my $copy_ok = $helpers->copy_file($source_path, $destination_path);
    return $copy_ok
}

1;

