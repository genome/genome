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
        user => {
            is => 'UR::Object',
            doc => 'To generate a SoftwareResult::User',
        },
        label => {
            is => 'Text',
            doc => 'Add the user with this label',
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
    my $class = $self->class;
    die "Abstract method 'input_hash' must be defined in class $class";
}

# Do whatever you want to after a result was created or looked up.
# Return undef on failure.
sub post_get_or_create {
    my $self = shift;
    return 1;
}

sub shortcut {
    my $self = shift;

    return $self->_protected_fetch('shortcut', 'get_with_lock', 'shortcut');
}

sub execute {
    my $self = shift;

    return $self->_protected_fetch('execute', 'get_or_create', 'created');
}

sub _protected_fetch {
    my ($self, $name, $method, $label) = validate_pos(@_, {type => OBJECT},
        {type => SCALAR}, {type => SCALAR}, {type => SCALAR});

    my $error;
    my $rv = try {
        $self->_fetch_result($method, $label);
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
    my ($self, $method, $user_label) = validate_pos(@_, {type => OBJECT},
        {type => SCALAR}, {type => SCALAR});

    $self->debug_message("Attempting to %s a %s with arguments %s",
        $method, $self->result_class, pp({$self->input_hash}));
    my $result = $self->result_class->$method($self->input_hash);

    if ($result) {
        $self->debug_message("%s returned result (%s)", $method, $result->id);
        $self->output_result($result);
        $self->create_software_result_user($user_label);
        if (defined $self->label) {
            $self->create_software_result_user($self->label);
        }
        return $self->post_get_or_create;
    } else {
        $self->debug_message("Failed to %s result.", $method);
        return 0;
    }

}

sub create_software_result_user {
    my ($self, $label) = @_;

    if ($self->user) {
        $self->debug_message(
            "Making %s(%s) a user of SoftwareResult(%s) with label '%s'",
            $self->user->class, $self->user->id, $self->output_result->id, $label,
        );
        $self->output_result->add_user(user => $self->user, label => $label);
    } else {
        $self->debug_message("Not making anything a user of " .
            "SoftwareResult(%s) because no 'user' was provided as an input",
            $self->output_result->id
        );
    }
    return;
}

1;
