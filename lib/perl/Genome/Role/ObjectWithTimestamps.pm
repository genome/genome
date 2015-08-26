package Genome::Role::ObjectWithTimestamps;

use strict;
use warnings;

use UR::Role;

role Genome::Role::ObjectWithTimestamps {
    has => [
        created_at => { is => 'Timestamp' },
        updated_at => { is => 'Timestamp' },
    ],
};

sub __import__ {
    my($role_name, $class_meta) = @_;

    my $class_name = $class_meta->class_name;
    UR::Observer->register_callback(
        subject_class_name => $class_name,
        callback => \&is_updated,
    );

    UR::Observer->register_callback(
        subject_class_name => $class_name,
        aspect => 'create',
        callback => \&is_created,
    );
}

sub is_created {
    my $self = shift;
    unless ($self->created_at) {
        $self->created_at(UR::Context->current->now);
    }
    return 1;
}

sub is_updated {
    my ($self, $aspect) = @_;
    return unless ref($self);  # only interested in notifications to instances
    if (ref($self) && $aspect ne 'commit' && $aspect ne 'load' && $aspect ne 'updated_at') {
        $self->updated_at(UR::Context->current->now);
    }
}

1;
