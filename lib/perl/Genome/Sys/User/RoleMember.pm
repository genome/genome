package Genome::Sys::User::RoleMember;

use strict;
use warnings;

use Genome;
use Carp;

class Genome::Sys::User::RoleMember {
    table_name => 'subject.role_member',
    id_by => [
        user_email => {
            is => 'Text',
            len => 255,
        },
        role_id => {
            is => 'Text',
            len => 32,
        },
    ],
    has => [
        user => {
            is => 'Genome::Sys::User',
            id_by => 'user_email',
            constraint_name => 'GSURM_UE_FK',
        },
        role => {
            is => 'Genome::Sys::User::Role',
            id_by => 'role_id',
            constraint_name => 'GSURM_RI_FK',
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

Genome::Sys::User::RoleMember->add_observer(
    callback => \&_change_callback,
);
sub _change_callback {
    my ($self, $signal) = @_;
    return 1 unless grep { $signal eq $_ } qw/ precommit create delete /;
    unless (Genome::Sys->current_user_is_admin) {
        Carp::confess 'Only admins can modify role/member bridges!';
    }
    return 1;
}

sub create {
    my $class = shift;
    my $bool_expr = UR::BoolExpr->resolve_normalized($class, @_);

    my $email = $bool_expr->value_for('user_email');
    unless ($email) {
        Carp::confess "Create method for $class not given user id, cannot resolve user object!";
    }
    my $user = Genome::Sys::User->get(email => $email);
    unless ($user) {
        Carp::confess "Could not resolve user from email $email!";
    }

    my $role_id = $bool_expr->value_for('role_id');
    unless ($role_id) {
        Carp::confess "Create method for $class not given role id, cannot resolve role object!";
    }
    my $role = Genome::Sys::User::Role->get($role_id);
    unless ($role) {
        Carp::confess "Could not resolve role from id $role_id!";
    }

    if ($user->has_role($role)) {
        Carp::confess "User " . $user->__display_name__ . " already has role " . $role->__display_name__ . "!";
    }

    return $class->SUPER::create(@_);
}

1;

