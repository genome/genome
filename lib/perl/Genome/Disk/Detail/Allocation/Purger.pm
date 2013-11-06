package Genome::Disk::Detail::Allocation::Purger;

use strict;
use warnings;

use Genome;

use Carp qw(confess);
use File::Copy::Recursive qw(dircopy dirmove);

class Genome::Disk::Detail::Allocation::Purger {
    is => 'Genome::Disk::Detail::StrictObject',

    has => [
        allocation_id => {
            is => 'Text',
            len => 64,
        },
        reason => {
            is => 'Text',
        },
    ],

};


sub purge {
    my $self = shift;

    my $allocation_object = Genome::Disk::Allocation->get($self->allocation_id);
    die sprintf("No allocation found for id: %s",
        $self->allocation_id) unless $allocation_object;

    $allocation_object->_create_file_summaries();
    my $trash_folder = $allocation_object->_get_trash_folder();

    unless($ENV{UR_DBI_NO_COMMIT}) {
        my $destination_directory = Genome::Sys->create_directory(
            File::Spec->join($trash_folder, $allocation_object->id));
        dirmove($allocation_object->absolute_path, $destination_directory);
    }

    Genome::Timeline::Event::Allocation->purged(
        $self->reason,
        $allocation_object,
    );

    $allocation_object->status('purged');
    $allocation_object->kilobytes_requested(0);
    $allocation_object->kilobytes_used(0);

    return 1;
}


1;
