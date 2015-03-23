package Genome::Ptero::Wrapper;

use strict;
use warnings FATAL => qw(all);
use Genome;
use Ptero;
use UR;
use Data::Dump qw();

class Genome::Ptero::Wrapper {
    is => 'Command::V2',
    has => [
        'command_class' => {
            is => 'Text',
            doc => 'class name of the Genome Command to execute',
        },
        'method' => {
            is => 'Text',
            valid_values => ['execute', 'shortcut'],
            doc => 'method to call on the Genome Command object',
        },
    ],

    doc => 'Ptero execution wrapper for Genome Command objects',
};


sub execute {
    my $self = shift;

    my $execution = Ptero::Workflow::Execution->new();

    $self->_setup_logging($execution->metadata);
    $self->_log_execution_information($execution);

    my $command = $self->_instantiate_command($execution->inputs);

    $self->_run_command($command);
    $execution->set_outputs(
        _get_command_outputs($command, $self->command_class));

    return 1;
}

sub _setup_logging {
    my ($self, $execution_metadata) = @_;
    # Replace stdout/stderr with files
}

sub _log_execution_information {
    my ($self, $execution) = @_;
    # - Dump entire execution response
    # - Write a pastable command that sets the PTERO_WORKFLOW_EXECUTION_URL and
    #   runs this wrapper (could be useful for debugging).
}

sub _instantiate_command {
    my ($self, $inputs) = @_;

    my $pkg = $self->command_class;
    eval "use $pkg";
    my $cmd = eval {$pkg->create(%$inputs)};
    my $error = $@;
    if ($error) {
        Carp::confess sprintf(
            "Failed to intantiate class (%s) with inputs (%s): %s",
            $pkg, Data::Dump::pp($inputs), $error);
    }
    return $cmd;
}

sub _run_command {
    my ($self, $command) = @_;

    my $method = $self->method;
    my $ret = eval { $command->$method() };
    my $error = $@;
    if ($error) {
        Carp::confess sprintf(
            "Crashed in %s for command %s: %s",
            $self->method, $self->command_class, $error,
        );
    }
    unless ($ret) {
        Carp::confess sprintf("Failed to %s for command %s.",
            $self->method, $self->command_class,
        );
    }

    _commit();

    $self->status_message("Succeeded to %s command %s",
        $method, $self->command_class);
}

sub _commit {
    eval {
        UR::Context->commit();
    };

    my $error = $@;
    if ($error) {
        Carp::confess "Failed to commit: $error";
    }
}

sub _get_command_outputs {
    my ($cmd, $pkg) = @_;

    my %outputs;
    for my $prop (_output_names($pkg)) {
        my $prop_name = $prop->property_name;
        my $value = $prop->is_many ? [$cmd->$prop_name] : $cmd->$prop_name;
        $outputs{$prop_name} = $value;
    }

    $outputs{result} = 1 unless exists $outputs{result};

    return \%outputs;
}

sub _output_names {
    my $pkg = shift;

    return $pkg->__meta__->properties(is_output => 1);
}


1;
