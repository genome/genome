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

# Eventually SoftwareResults should have something like a status field,
# removing the need to set a test_name to invalidate a SoftwareResult
# backed by a purged allocation.
sub _update_owner_test_name {
    my $self       = shift;
    my $allocation = shift;
    my $event      = shift;

    my $owner = $allocation->owner;

    # Make sure the owner is a SoftwareResult before we set_test_name.
    return unless defined $owner;
    return unless $owner->isa('Genome::SoftwareResult');
    return if defined $owner->test_name;

    my $timestamp = $event->created_at;
    my $reason    = $event->reason;

    my $test_name = sprintf "Allocation ID %s purged on %s, %s",
        $allocation->id, $timestamp,
        (defined $reason ? "with reason '$reason'" : 'no reason specified');

    $owner->set_test_name($test_name);

    return;
}

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

    my $event = Genome::Timeline::Event::Allocation->purged(
        $self->reason,
        $allocation_object,
    );

    $self->_update_owner_test_name($allocation_object, $event);

    $allocation_object->status('purged');
    $allocation_object->kilobytes_requested(0);
    $allocation_object->kilobytes_used(0);

    return 1;
}


1;
