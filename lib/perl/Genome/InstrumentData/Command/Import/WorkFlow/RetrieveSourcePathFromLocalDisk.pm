package Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathFromLocalDisk;

use strict;
use warnings;

use Genome;

require File::Copy;

class Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathFromLocalDisk {
    is => 'Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePath',
};

sub _path_size {
    return -s $_[0];
}

sub _source_path_size {
    return _path_size($_[0]->source_path);
}

sub _retrieve_source_path {
    my $self = shift;

    my $copy_ok = File::Copy::copy($self->source_path, $self->destination_path);
    if ( not $copy_ok ) {
        $self->error_message('Copy from %s to %s failed!', $self->source_path, $self->destination_path);
        return;
    }

    return 1;
}

sub source_md5 {
    my $self = shift;

    my $md5_path = $self->helpers->md5_path_for($self->source_path);
    my $md5_path_exists = _path_size($md5_path);
    return if not $md5_path_exists;

    return $self->helpers->load_md5($md5_path);
}

1;

