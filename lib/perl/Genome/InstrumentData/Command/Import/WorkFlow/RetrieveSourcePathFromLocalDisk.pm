package Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathFromLocalDisk;

use strict;
use warnings;

use Genome;

require File::Basename;

class Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathFromLocalDisk {
    is => 'Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePath',
};

sub _retrieve_path {
    my ($self, $from, $to) = @_;
    $self->debug_message('Retrieve path via copy...');

    $self->debug_message("From: $from");
    $self->debug_message("To: $to");

    my $move_ok = File::Copy::copy($from, $to);
    if ( not $move_ok ) {
        $self->error_message('Copy failed!');
        return;
    }

    my $from_sz = -s $from;
    $self->debug_message("From size: $from_sz");
    my $to_sz = -s $to;
    $self->debug_message("To size: $to_sz");
    if ( not $to_sz or $to_sz != $from_sz ) {
        $self->error_message("Copy succeeded, but destination size is different from original! $to_sz vs $from_sz");
        return;
    }

    $self->debug_message('Retrieve path via copy...done');
    return 1;
}

1;

