package Genome::Command::DelegatesToResult;

use strict;
use warnings FATAL => 'all';
use Genome;
use Data::Dump qw(pp);
use Params::Validate qw(validate_pos :types);
use Try::Tiny qw(try catch);

class Genome::Command::DelegatesToResult {
    is => ['Command::V2'],
    has_output => [
        output_result => {
            is => 'Genome::SoftwareResult',
        },
    ],
    has_optional_input => [
        sponsor => {
            is => 'UR::Object',
            doc => 'The responsible party for this result (defaults to the current user)',
        },
        requestor => {
            is => 'UR::Object',
            doc => 'The process that wants the result (must be provided if not derivable)',
        },
        user => {
            is => 'UR::Object',
            doc => 'An additional user of the result',
        },
        label => {
            is => 'Text',
            doc => 'Label for the additional user (if any) of the result',
        }
    ],
};

sub result_class {
    my $self = shift;
    my $class = $self->class;
    die "Abstract method 'result_class' must be defined in class $class";
}

sub input_hash {
    my $self = shift;

    my %inputs = ();
    for my $property ($self->result_class->__meta__->properties()) {
        next unless $property->{is_param} or $property->{is_input};

        my $name = $property->property_name;
        next unless $self->can($name);

        if($property->is_many) {
            $inputs{$name} = [$self->$name];
        } else {
            $inputs{$name} = $self->$name;
        }
    }

    return %inputs;
}

sub _input_hash {
    my $self = shift;

    my $result_users = {
        sponsor   => $self->sponsor // Genome::Sys->current_user,
        requestor => $self->requestor,
    };
    if($self->label and $self->user) {
        $result_users->{$self->label} = $self->user;
    }

    return (
        test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME},
        $self->input_hash,
        users => $result_users,
    );
}

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__;

    for my $pair (['label', 'user'], ['user', 'label']) {
        my ($first_property, $second_property) = @$pair;
        if(defined $self->$first_property and not defined $self->$second_property) {
            push @errors, UR::Object::Tag->create(
                type => 'error',
                properties => [$second_property],
                desc => sprintf(
                    q('%s' is required when '%s' is specified.),
                    $second_property, $first_property
                ),
            );
        }
    }

    unless(defined $self->requestor) {
        push @errors, UR::Object::Tag->create(
            type => 'error',
            properties => 'requestor',
            desc => 'No value specified for required property',
        );
    }

    return @errors;
}

# Do whatever you want to after a result was created or looked up.
# Return undef on failure.
sub post_get_or_create {
    my $self = shift;
    return 1;
}

sub shortcut {
    my $self = shift;

    return $self->_protected_fetch('shortcut', 'get_with_lock');
}

sub execute {
    my $self = shift;

    return $self->_protected_fetch('execute', 'get_or_create');
}

sub _protected_fetch {
    my ($self, $name, $method) = validate_pos(@_, {type => OBJECT},
        {type => SCALAR}, {type => SCALAR});

    my $error;
    my $rv = try {
        $self->_fetch_result($method);
    } catch {
        $error = $_ || 'unknown error';
    };

    if ($error) {
        $self->error_message("Exception in %s: %s", $name, $error);
        return;
    } else {
        return $rv;
    }
}

sub _fetch_result {
    my ($self, $method) = validate_pos(@_, {type => OBJECT}, {type => SCALAR});

    $self->debug_message("Attempting to %s a %s with arguments %s",
        $method, $self->result_class, pp({$self->_input_hash}));
    my $result = $self->result_class->$method($self->_input_hash);

    if ($result) {
        $self->debug_message("%s returned result (%s)", $method, $result->id);
        $self->delete_transients($result);
        $self->output_result($result);
        return $self->post_get_or_create;
    } else {
        $self->debug_message("Failed to %s result.", $method);
        return 0;
    }

}

sub delete_transients {
    my ($self, $result) = @_;
    for my $transient_property ($result->__meta__->properties(is_transient => 1)) {
        next if $transient_property->property_name eq '_lock_name';

        my $value = undef;
        if (defined($transient_property->default_value)) {
            $value = $transient_property->default_value;
        }

        my $name = $transient_property->property_name;
        $result->$name($value);
    }
    return;
}

1;
