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
    ],

    has_optional => [
        id => {
            is => 'Text',
            len => 64,
        },

        group_subdirectory => {
            is => 'Text',
            len => 255,
        },

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

sub get_id {
    my $self = shift;
    if ($self->id) {
        return $self->id;
    } else {
        # TODO autogenerate_new_object_id should technically receive a BoolExpr
        return Genome::Disk::Allocation->__meta__->autogenerate_new_object_id;
    }
}

sub sanitize {
}

# TODO This needs to be removed, site-specific
our @APIPE_DISK_GROUPS = qw/
    info_apipe
    info_apipe_ref
    info_alignments
    info_genome_models
    research
    systems_benchmarking
/;
sub validate {
    my $self = shift;

    $self->validate_owner_class_name;
    $self->validate_kilobytes_requested;
    $self->validate_disk_group_name;
}

sub validate_owner_class_name {
    my $self = shift;

    unless ($self->owner_class_name->__meta__) {
        confess sprintf("Could not find meta information for owner class %s, "
            . "make sure this class exists!", $self->owner_class_name);
    }
}

sub validate_kilobytes_requested {
    my $self = shift;

    unless ($self->kilobytes_requested >= 0) {
        confess sprintf('Kilobytes requested is negative (%s)!',
            $self->kilobytes_requested);
    }
}

sub validate_disk_group_name {
    my $self = shift;

    unless (grep { $self->disk_group_name eq $_ } @APIPE_DISK_GROUPS) {
        confess "Can only allocate disk in apipe disk groups, "
            . "not %s. Apipe groups are: %s",
            $self->disk_group_name, join(", ", @APIPE_DISK_GROUPS);
    }
}


1;
