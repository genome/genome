package Genome::Disk::Detail::Allocation::CreationParameters;

use strict;
use warnings;

use Genome;
use UR;

use Carp qw(confess);

class Genome::Disk::Detail::Allocation::CreationParameters {
    has => [
        kilobytes_requested => {
            is => 'Number',
            len => 20,
        },

        owner_class_name => {
            is => 'Text',
            len => 255,
        },
        owner_id => {
            is => 'Text',
            len => 255,
        },

        allocation_path => {
            is => 'Text',
            len => 4000,
        },
        disk_group_name => {
            is => 'Text',
            len => 40,
        },

        group_subdirectory => {
            is => 'Text',
            len => 255,
        },
    ],

    has_optional => [
        mount_path => {
            is => 'Text',
            len => 255,
        },
        exclude_mount_path => {
            is => 'Text',
            len => 255,
        },

        archive_after_time => {
            is => 'DateTime',
            len => 11,
        },
        kilobytes_used => {
            is => 'Number',
            len => 20,
        },
    ],
};


sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);

    if ($self->__errors__) {
        my @messages;
        for my $error_tag ($self->__errors__) {
            push @messages, $error_tag->__display_name__ . "\n";
        }
        confess sprintf(
            'Could not create params object:\n%s', join('', @messages));
    }

    $self->sanitize;
    $self->validate;

    return $self;
}

sub sanitize {
}

sub validate {
}


1;
