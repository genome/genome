package Genome::Command::DelegatesToResult;

use strict;
use warnings FATAL => 'all';
use Genome;
use Data::Dump qw(pp);

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
sub post_get_or_create {
    my $self = shift;
    return;
}

sub shortcut {
    my $self = shift;

    $self->debug_message("Attempting to get a %s with arguments %s",
        $self->result_class, pp({$self->input_hash}));
    my $result = $self->result_class->get_with_lock($self->input_hash);
    if ($result) {
        $self->debug_message("Found existing result (%s)", $result->id);
        $self->output_result($result);
        $self->create_software_result_user('shortcut');
        $self->post_get_or_create;
        return 1;
    } else {
        $self->debug_message("Found no existing result.");
        return 0;
    }
}

sub execute {
    my $self = shift;

    $self->debug_message("Attempting to get or create a %s with arguments %s",
        $self->result_class, pp({$self->input_hash}));
    my $result = $self->result_class->get_or_create($self->input_hash);

    if ($result) {
        $self->debug_message("Got or created result (%s)", $result->id);
        $self->output_result($result);
        $self->create_software_result_user('created');
        $self->post_get_or_create;
        return 1;
    } else {
        $self->debug_message("Failed to get or create result.");
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
        return;
    } else {
        $self->debug_message("Not making anything a user of " .
            "SoftwareResult(%s) because no 'user' was provided as an input",
            $self->output_result->id
        );
        return;
    }
}

1;
