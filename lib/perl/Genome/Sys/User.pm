package Genome::Sys::User;

use strict;
use warnings;

use Genome;
use Carp;

class Genome::Sys::User {
    is => 'Genome::Searchable',
    table_name => 'subject.user',
    id_by => [
        email => {
            is => 'Text',
            len => 255,
            doc => 'Email of the user, must be unique',
        },
    ],
    has_optional => [
        name => {
            is => 'Text',
            len => 64,
            doc => 'Full name of the user (eg, Ronald McDonald)',
        },
        username => {
            is => 'Text',
            len => 64,
            doc => 'System user name of the user (eg, rmcdonald)',
        },
    ],
    has_many_optional => [
        project_parts => {
            is => 'Genome::ProjectPart',
            reverse_as => 'entity',
            is_mutable => 1,
        },
        projects => {
            is => 'Genome::Project',
            via => 'project_parts',
            to => 'project',
            is_mutable => 1,
        },
        project_names => {
            is => 'Text',
            via => 'projects',
            to => 'name',
        },
        user_role_bridges => {
            is => 'Genome::Sys::User::RoleMember',
            reverse_as => 'user',
        },
        user_roles => {
            is => 'Genome::Sys::User::Role',
            via => 'user_role_bridges',
            to => 'role',
            is_mutable => 1,
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

sub _resolve_param_value_from_text_by_name_or_id {
    my $class = shift;
    my $param_arg = shift;

    #First try default behaviour of looking up by name or id
    my @results = Genome::Command::Base->_resolve_param_value_from_text_by_name_or_id($class, $param_arg);

    unless (@results) {
        @results = $class->get(username => $param_arg)
    }

    unless (@results) {
        @results = $class->get(name => $param_arg)
    }

    return @results;
}

Genome::Sys::User->add_observer(
    callback => \&_change_callback,
);
sub _change_callback {
    my ($self, $signal) = @_;
    return 1 unless grep { $signal eq $_ } qw/ create delete precommit /;
    unless (Genome::Sys->current_user_is_admin) {
        Carp::confess 'Only admins can modify users!';
    }
    return 1;
}

sub create {
    my $class = shift;
    my $bool_expr = UR::BoolExpr->resolve_normalized($class, @_);

    my $email = $bool_expr->value_for('email');
    unless ($email) {
        Carp::confess "Cannot create a user without an email!";
    }
    my @users = Genome::Sys::User->get(email => $email);
    if (@users) {
        if (@users == 1) {
            Carp::confess "Another user already exists with email address $email, cannot create another!";
        }
        else {
            Carp::confess "Somehow there are many users with email address $email... Cannot create another, please contact informatics!";
        }
    }

    return $class->SUPER::create(@_);
}

sub delete {
    my $self = shift;
    for my $bridge ($self->user_role_bridges) {
        my $role_name = $bridge->role->name;
        my $user_name = $self->name;
        my $rv = $bridge->delete;
        unless ($rv) {
            Carp::confess "Could not delete bridge object between user $user_name and role $role_name";
        }
    }
    return $self->SUPER::delete(@_);
}

sub fix_params_and_get {
    my ($class, @p) = @_;
    my %p;
    if (scalar(@p) == 1) {
        my $key = $p[0];
        $p{'email'} = $key;
    }
    else {
        %p = @p;
    }

    if (defined($p{'email'})
        && $p{'email'} !~ /\@/) {
        my $old = $p{'email'};
        my $new = join('@',$p{'email'},Genome::Config::domain());
        warn "Trying to get() for '$old' - assuming you meant '$new'";
        $p{'email'} = $new;
    }

    return $class->SUPER::get(%p);
}

sub has_role_by_name {
    my ($self, $role_name) = @_;
    my $role = Genome::Sys::User::Role->get(name => $role_name);
    return $self->has_role($role);
}

sub has_role {
    my ($self, $role) = @_;
    return 0 unless $role;
    my $bridge = Genome::Sys::User::RoleMember->get(
        user => $self,
        role => $role,
    );
    return defined $bridge;
}

1;
