package Genome::Disk::Command::Allocation::Import;

use strict;
use warnings;

use Genome;
use Cwd 'realpath';

class Genome::Disk::Command::Allocation::Import {
    is => 'Command::V2',
    has => [
        path => {
            is => 'DirectoryPath',
            doc => 'Path at which data exists, allocation will be made to point to this location',
        },
        owner_class_name => {
            is => 'Text',
            default_value => 'Genome::Sys::User',
            doc => 'Class name of entity that owns the allocations (eg, Genome::Sys::User for users, or Genome::Model::Build::* for builds.',
        },
        owner_id => {
            is => 'Text',
            default_value => Genome::Sys::User->owner_id,
            doc => 'The ID used to retrieve the owner (in conjunction with owner_class_name), ' .
                   'e.g. username@example.com for Genome::Sys::User or build_id for Genome::Model::Build.',
        },
    ],
    has_optional => [
        volume_prefix => {
            is => 'Text',
            default_value => '/gscmnt',
            doc => 'Volume prefix to expect at start of the path',
        },
        allocation => {
            is => 'Genome::Disk::Allocation',
            is_transient => 1,
            doc => 'Transient argument, not settable via command line. Useful for programmatic access to created allocation',
        },
        temp_allocation_path => {
            is => 'DirectoryPath',
            is_transient => 1,
            doc => 'Transient argument, not settable via command line. Stores temporary location that allocation is original set to',
        },
    ],
    doc => 'Creates an allocation for data that already exists on disk',
};

sub help_brief {
    return 'Creates an allocation for data that already exists on disk';
}

sub help_synopsis {
    return help_brief();
}

sub help_detail {
    return <<EOS
This tool creates an allocation for data that already exists. The 'typical'
allocation create process assumes that data does not already exist at the
location specified and will fail if this is not true. This command
works around this issue.
EOS
}

sub execute {
    my $self = shift;

    my $kb = $self->_validate_and_get_size_of_path($self->path);
    my ($mount_path, $group_subdir, $allocation_path) = $self->_parse_path($self->path, $self->volume_prefix);
    my $volume = $self->_get_and_validate_volume($mount_path);
    my $group = $self->_get_and_validate_group($group_subdir, $volume->groups);
    $self->_verify_allocation_path($allocation_path);

    $self->status_message("Parsed and verified target path " . $self->path . ", creating allocation.");

    my %params = (
        disk_group_name => $group->disk_group_name,
        allocation_path => sprintf("%s-temp_allocation_path_for_existing_data", $allocation_path),
        kilobytes_requested => $kb,
        owner_class_name => $self->owner_class_name,
        owner_id => $self->owner_id,
        mount_path => $mount_path,
    );
    my $allocation = Genome::Disk::Allocation->create(%params);
    unless ($allocation) {
        require Data::Dumper;
        die "Could not create allocation with these params:\n" . Data::Dumper::Dumper(\%params);
    }
    my $temp_location = $allocation->absolute_path;
    $self->temp_allocation_path($temp_location);

    $self->status_message("Created allocation " . $allocation->id . ", moving to final location");

    $allocation->allocation_path($allocation_path);
    unless ($allocation->absolute_path eq $self->path) {
        die "Somehow, absolute path of new allocation does not match expected value " . $self->path;
    }

    unless (rmdir($temp_location)) {
        $self->warning_message("Could not remove temporary path $temp_location");
    }

    $self->allocation($allocation);
    return 1;
}

sub _verify_allocation_path {
    my ($class, $allocation_path) = @_;
    # I don't know why this method is private, I don't think it should be
    if (my $parent_allocation = Genome::Disk::Allocation->get_parent_allocation($allocation_path)) {
        die "Parent allocation (" . $parent_allocation->id . ") found for $allocation_path!";
    }

    # Ditto... don't think this should be private
    if (Genome::Disk::Allocation->get_child_allocations($allocation_path)) {
        die "Child allocation(s) found for $allocation_path!";
    }

    return 1;
}

sub _parse_path {
    my ($class, $target_path, $volume_prefix) = @_;
    die "_parse_path method must be given a target path!" unless $target_path;
    die "_parse_path method must be given a volume prefix!" unless $volume_prefix;
    my ($mount_path, $group_subdir, $allocation_path) =
        $target_path =~  /^($volume_prefix\/\w+)\/(\w+)\/(.+)$/;
    unless ($mount_path and $group_subdir and $allocation_path) {
        die "Could not determine mount path, group subdirectory, or allocation path from given path!";
    }
    return ($mount_path, $group_subdir, $allocation_path);
}

sub _get_and_validate_volume {
    my ($class, $mount_path) = @_;
    my $volume = Genome::Disk::Volume->get(
        mount_path => $mount_path,
        disk_status => 'active',
        can_allocate => 1,
    );
    unless ($volume) {
        die "Found no allocatable and active volume with mount path $mount_path";
    }
    return $volume;
}

sub _get_and_validate_group {
    my ($class, $group_subdir, @groups) = @_;
    my $group;
    for my $candidate_group (@groups) {
        if ($candidate_group->subdirectory eq $group_subdir) {
            $group = $candidate_group;
            last;
        }
    }
    unless ($group) {
        die "No groups found with subdirectory that matches $group_subdir";
    }
    return $group;
}

sub _validate_and_get_size_of_path {
    my ($class, $target_path) = @_;
    if (-l $target_path) {
        die "Path is a link, you should try allocating " . realpath($target_path)
            . " instead of " . $target_path;
    }
    if (not -d $target_path) {
        die "Path does not exist or is not a directory: " . $target_path;
    }
    if (realpath($target_path) ne $target_path) {
        die "Target path and absolute path differ; you should try allocating "
            . realpath($target_path) . " instead of " . $target_path;
    }
    my $kb = Genome::Sys->disk_usage_for_path($target_path);
    unless (defined $kb) {
        die "Could not determine size of " . $target_path;
    }
    return $kb;
}

1;
