package Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePaths;

use strict;
use warnings;

use Genome;

require File::Basename;

class Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePaths { 
    is => 'Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathBase',
    has_input => [
        source_paths => {
            is => 'Text',
            is_many => 1,
            doc => 'Source paths of sequences to get.',
        },
    ],
    has_output => [
        destination_paths => {
            is => 'Text',
            is_many => 1,
            calculate => q| return map { $self->destination_path_for_source_path($_) } $self->source_paths; |,
            doc => 'Final destination paths.',
        }, 
    ],
};

sub execute {
    my $self = shift;

    for my $source_path ( $self->source_paths ) {
        my $retrieve_source_path = $self->_retrieve_source_path($source_path);
        return if not $retrieve_source_path;

        my $retrieve_source_md5 = $self->_retrieve_source_md5($source_path);
        return if not $retrieve_source_md5;
    }

    return 1;
}

1;

