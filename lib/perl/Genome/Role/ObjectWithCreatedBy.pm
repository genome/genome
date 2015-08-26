package Genome::Role::ObjectWithCreatedBy;

use strict;
use warnings;

use UR::Role;

role Genome::Role::ObjectWithCreatedBy {
    has => [
        created_by => { is => 'Text' },
    ],
};

sub __import__ {
    my($role_name, $class_meta) = @_;

    my $class_name = $class_meta->class_name;
    UR::Observer->register_callback(
        subject_class_name => $class_name,
        aspect => 'create',
        callback => \&_populate_created_by
    );
}

sub _populate_created_by {
    my $self = shift;
    unless ($self->created_by) {
        $self->created_by(Genome::Sys->username);
    }
}

1;
