package Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathFromLocalDisk;

use strict;
use warnings;

use Genome;

require File::Copy;

class Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathFromLocalDisk {
    is => 'Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePath',
};

sub _retrieve_source_path {
    my $self = shift;

    my $copy_ok = File::Copy::copy($self->source_path, $self->destination_path);
    if ( not $copy_ok ) {
        $self->error_message('Copy from %s to %s failed!', $self->source_path, $self->destination_path);
        return;
    }

    return 1;
}

sub _load_source_md5 {
    my $self = shift;
    return if not $self->source_file->md5_path_size;
    return Genome::InstrumentData::Command::Import::WorkFlow::Helpers->load_md5($self->source_file->md5_path);
}

1;

