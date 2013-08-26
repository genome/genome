package Genome::Sys::User::Role;

use strict;
use warnings;

use Genome;
use Carp;

class Genome::Sys::User::Role {
    table_name => 'subject.role',
    id_by => [
        id => {
            is => 'Text',
            len => 32,
        },
    ],
    has => [
        name => {
            is => 'Text',
            len => 64,
        },
    ],
    has_many_optional => [
        user_bridges => {
            is => 'Genome::Sys::User::RoleMember',
            reverse_as => 'role',
        },
        users => {
            is => 'Genome::Sys::User',
            via => 'user_bridges',
            to => 'user',
            is_mutable => 1,
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    id_generator => '-uuid',
};

Genome::Sys::User::Role->add_observer(
    callback => \&_change_callback,
);

sub __display_name__ {
    my $self = shift;
    return $self->name;
}

sub _change_callback {
    my ($self, $signal) = @_;
    return 1 unless grep { $signal eq $_ } qw/ create delete precommit /;
    unless (Genome::Sys->current_user_is_admin) {
        Carp::confess "Only admins can change role names!";
    }
    return 1;
}

sub create {
    my $class = shift;
    my $bool_expr = UR::BoolExpr->resolve_normalized($class, @_);

    my $name = $bool_expr->value_for('name');
    unless ($name) {
        Carp::confess "Cannot create a role without a name!";
    }
    my @roles = Genome::Sys::User::Role->get(
        name => $name,
    );
    if (@roles) {
        if (@roles == 1) {
            Carp::confess "Another role with name $name already exists, cannot create another!";
        }
        else {
            Carp::confess "Somehow there are " . scalar @roles . " roles with name $name. Cannot create another role, please contact informatics about this...";
        }
    }

    return $class->SUPER::create(@_);
}

sub delete {
    my $self = shift;
    my @users = $self->users;
    if (@users) {
        Carp::confess "Cannot delete user role " . $self->name .
            ", the following users use it: " . join(', ', map { $_->name } @users);
    }

    return $self->SUPER::delete(@_);
}

1;

