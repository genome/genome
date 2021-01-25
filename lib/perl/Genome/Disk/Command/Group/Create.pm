package Genome::Disk::Command::Group::Create;

use strict;
use warnings;

use Genome;
use Genome::Utility::Text qw();

class Genome::Disk::Command::Group::Create {
    is => 'Command::V2',
    doc => 'Create a new disk group',
    has_input => [
        name => {
            is => 'Text',
            doc => 'name of the new disk group',
        },
        subdirectory => {
            is => 'Text',
            doc => 'the directory within a volume under which all allocations will be created',
        },
        unix_group => {
            is => 'Text',
            doc => 'the group who owns this disk group',
        },
    ],
};

sub execute {
    my $self = shift;

    my $name = $self->name;
    my $subdirectory = $self->subdirectory;
    my $unix_group = $self->unix_group;

    my $existing = Genome::Disk::Group->get(name => $name);
    if ($existing) {
        $self->fatal_message('A group named %s already exists.', $name);
    }

    my $gid = Genome::Sys::gidgrnam($unix_group);

    my $sanitized_subdirectory = Genome::Utility::Text::sanitize_string_for_filesystem($subdirectory);
    unless ($sanitized_subdirectory eq $subdirectory) {
        $self->fatal_message('Proposed subdirectory %s contains invalid characters.', $subdirectory);
    }

    my $new_group = Genome::Disk::Group->create(
        name => $name,
        subdirectory => $subdirectory,
        setgid => 0,
        unix_uid => 0,
        unix_gid => $gid,
        permissions => 770,
    );

    unless ($new_group) {
        $self->fatal_message('Failed to create group.');
    }
    $self->status_message('New group %s created with ID %s.', $name, $new_group->id);

    return $new_group->id;
}

1;
