package Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePath;

use strict;
use warnings;

use Genome;

require File::Basename;

class Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePath { 
    is => 'Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathBase',
    has_input => [
        source_path => {
            is => 'Text',
            doc => 'Source path of sequences to get.',
        },
    ],
    has_output => [
        destination_path => {
            calculate => q| return $self->destination_path_for_source_path($self->source_path); |,
            doc => 'Final destination path.',
        }, 
    ],
};

sub execute {
    my $self = shift;

    my $source_path = $self->source_path;

    my $retrieve_source_path = $self->_retrieve_source_path($source_path);
    return if not $retrieve_source_path;

    my $retrieve_source_md5 = $self->_retrieve_source_md5($source_path);
    return if not $retrieve_source_md5;

    return 1;
}

1;

