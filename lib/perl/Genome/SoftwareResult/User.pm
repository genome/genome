package Genome::SoftwareResult::User;

use strict;
use warnings;
use Genome;
use List::MoreUtils qw(any);
use Params::Validate qw(:types);
use Carp qw();

use Genome::Utility::Text;

class Genome::SoftwareResult::User {
    table_name => 'result.user',
    id_by => [
        id => { is => 'Text', len => 32 },
    ],
    has => [
        software_result => {
            is => 'Genome::SoftwareResult',
            id_by => 'software_result_id',
            constraint_name => 'SRU_SR_FK',
        },
        user_class => {
            is => 'UR::Object::Type',
            id_by => 'user_class_name',
        },
        user_id => { is => 'Text', len => 256 },
        user => {
            is => 'UR::Object',
            id_by => 'user_id',
            id_class_by => 'user_class_name',
        },
        active => {
            is => 'Boolean',
            len => 1,
            default_value => 1,
            doc => 'Results actively being used should not be deleted',
        },
        label => { is => 'Text' },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    id_generator => '-uuid',
    doc => 'links a software result to other entities which depend on it',
};

sub with_registered_users {
    my $class = shift;
    my %params = Params::Validate::validate_with(
        params => \@_,
        spec   => {
            users => {
                type      => HASHREF,
                optional  => 0,
                callbacks => {
                    'must contain sponsor and requestor' => \&_validate_user_hash
                },
            },
            callback => {
                type     => CODEREF,
                optional => 0,
            }
        },
        allow_extra => 1,
    );

    my $user_hash = delete $params{users};
    my $sr_callback = delete $params{callback};

    my ($software_result, $newly_created) = $sr_callback->(%params);
    return unless $software_result;

    $class->_register_users($software_result, $user_hash, $newly_created);

    return $software_result;
}

sub _validate_user_hash {
    my $user_hash = shift;

    for my $type (qw( sponsor requestor )) {
        my $obj = $user_hash->{$type};
        return 0 unless($obj);
        return 0 unless(UNIVERSAL::isa($obj, _role_for_type($type)));
    }

    return 1;
}

sub _register_users {
    my $class = shift;
    my $software_result = shift;
    my $user_hash = shift;
    my $newly_created = shift;

    my %user_hash = %$user_hash;

    my $requestor = delete $user_hash{requestor};
    my $label = $newly_created ? 'created' : 'shortcut';
    $user_hash{$label} = $requestor;

    my @all_params;
    while(my ($label, $object) = each %user_hash) {
        my %params = (
            label           => $label,
            user            => $object,
            software_result => $software_result,
        );
        push @all_params, \%params;
    }

    my $observer;
    $observer = UR::Context->process->add_observer(
        aspect => 'precommit',
        callback => sub {
            if($observer) {
                $observer->delete;
                $observer = undef;
            } else {
                Carp::confess 'observer triggered multiple times!';
            }
            my @locks;
            for my $params (@all_params) {
                next if grep { $params->{$_}->isa('UR::DeletedRef') } qw(user software_result);
                push @locks, $class->_get_or_create_with_lock($params);
            }

            return unless @locks;

            my $unlocker;
            $unlocker = UR::Context->process->add_observer(
                aspect => 'commit',
                callback => sub {
                    if($unlocker) {
                        $unlocker->delete;
                        $unlocker = undef;
                    } else {
                        Carp::confess 'unlocker triggered multiple times!';
                    }

                    for my $lock (@locks) {
                        Genome::Sys->unlock_resource(resource_lock => $lock);
                    }
                }
            );
        }
    );
    unless($observer) {
        die 'Failed to create observer';
    }
}

sub _get_or_create_with_lock {
    my $class = shift;
    my $params = shift;

    my $existing = $class->get(%$params);
    return if $existing;

    my $resource = $class->_resolve_lock_name($params);
    my $lock = Genome::Sys->lock_resource(resource_lock => $resource, scope => 'site');
    $class->_get_or_create($params);

    return $lock;
}

sub _get_or_create {
    my $class = shift;
    my $params = shift;

    my $self = $class->load(%$params);
    unless($self) {
        $self = $class->create(%$params);
    }

    return $self;
}

sub _resolve_lock_name {
    my $class = shift;
    my $params = shift;

    return 'genome/software-result-user/' . Genome::Utility::Text::sanitize_string_for_filesystem(
        join('_', $params->{label}, $params->{user}->id, $params->{software_result}->id)
    );
}

sub _role_for_type {
    return sprintf(
        'Genome::SoftwareResult::%s',
        ucfirst(shift)
    );
}

sub user_hash_for_build {
    my $class = shift;
    my $build = shift;

    return {
        requestor => $build,
        sponsor   => $build->model->analysis_projects // Genome::Sys::User->get(username => $build->model->run_as),
    };
}

1;

