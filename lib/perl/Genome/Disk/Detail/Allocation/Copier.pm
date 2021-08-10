package Genome::Disk::Detail::Allocation::Copier;

use strict;
use warnings;

use Genome;

use Carp qw(confess);

class Genome::Disk::Detail::Allocation::Copier {
    is => 'Genome::Disk::Detail::StrictObject',

    has => [
        allocation_id => {
            is => 'Text',
            len => 64,
        },
        output_dir => {
            is => 'Text',
            len => 255,
        },
    ],
};

sub copy {
    my $self = shift;

    my $output_dir = $self->output_dir;
    Genome::Sys->create_directory($output_dir);

    my ($allocation_object, $allocation_lock) = Genome::Disk::Allocation->get_with_lock($self->allocation_id);

    my $unlocker = Scope::Guard->new(sub { $allocation_lock->unlock if $allocation_lock });

    my $original_path = $allocation_object->is_archived ? $allocation_object->archive_path : $allocation_object->absolute_path;

    # make shadow allocation
    my %creation_params = $self->_get_copy_shadow_params($allocation_object);

    # The shadow allocation is just a way of keeping track of our temporary
    # additional disk usage during the copy.
    my $shadow_allocation = Genome::Disk::Allocation->shadow_get_or_create(%creation_params);
    my %rsync_params = (
        source_directory => $original_path,
        target_directory => $output_dir,
    );
    my $shadow_absolute_path = $shadow_allocation->absolute_path;

    # copy files to output_dir
    my $copy_rv = eval {
        Genome::Sys->rsync_directory(%rsync_params);
    };
    unless ($copy_rv) {
        my $copy_error_message = $@;
        if (-d $shadow_absolute_path) {
            if (Genome::Sys->remove_directory_tree($shadow_absolute_path)) {
                $shadow_allocation->delete;
            }
        }
        confess(sprintf(
            "Could not copy allocation %s from %s to %s: %s",
            $allocation_object->id, $original_path,
            $shadow_absolute_path, $copy_error_message));
    }

    Genome::Timeline::Event::Allocation->copied(
        sprintf("copied from %s to %s", $original_path,
                $output_dir),
        $allocation_object,
    );

    Genome::Disk::Allocation::_commit_unless_testing();

    $shadow_allocation->delete;

}

sub _get_copy_shadow_path {
    my $allocation_path = shift;
    return sprintf("%s-copy_allocation_destination", $allocation_path)
}

sub _get_copy_shadow_params {
    my ($self, $allocation) = @_;

    my %creation_parameters = (
        disk_group_name => $allocation->disk_group_name,
        kilobytes_requested => $allocation->kilobytes_requested,
        owner_class_name => "UR::Value",
        owner_id => "shadow_allocation",
        allocation_path => _get_copy_shadow_path($allocation->allocation_path),
    );

    return %creation_parameters;
}


1;
