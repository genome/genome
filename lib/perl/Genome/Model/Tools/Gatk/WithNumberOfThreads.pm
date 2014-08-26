package Genome::Model::Tools::Gatk::WithNumberOfThreads;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Gatk::WithNumberOfThreads {
    is => 'UR::Object',
    is_abstract => 1,
    attributes_have => {
        is_input => { is => 'Boolean', is_optional => 1 },        
    },
    has_optional_input => {
        number_of_threads => {
            is => 'Number',
            doc => 'Controls the number of data threads sent to the processor',
        },
    },
    has_optional_transient => {
        shellcmd_exit_code => { is => 'Text', },
    },
};

sub _execute_command {
    my $self = shift;

    my $before_execute = $self->_before_execute;
    return if not $before_execute;

    my $run_gatk_command = $self->run_gatk_command;
    return if not $run_gatk_command;

    my $after_execute = $self->_after_execute;
    return if not $after_execute;

    return 1;
}

sub _before_execute { return 1; }
sub _after_execute { return 1; }

sub build_gatk_command {
    my $self = shift;
    my $command = $self->_build_gatk_command;
    return $command if not $self->number_of_threads;
    $command .= ' -nt '.$self->number_of_threads;
    return $command;
}

sub _build_gatk_command {
    Carp::confess('Please implement _build_gatk_command in sub class! '.$_[0]->class);
}

sub run_gatk_command {
    my $self = shift;

    my $command = $self->build_gatk_command;
    my $rv = eval{ Genome::Sys->shellcmd(cmd => $command); };
    if ( $rv ) {
        $self->shellcmd_exit_code(0);
        return 1;
    }

    my $error_message = "$@";
    $self->error_message($error_message);
    my $exit_code = $self->_resolve_exit_code_from_shellcmd_error_message($error_message);
    $self->shellcmd_exit_code($exit_code);
    $self->error_message("Failed to execute $command");

    return;
}

sub _resolve_exit_code_from_shellcmd_error_message {
    my ($self, $error_message) = @_;

    if ( $error_message =~ /^ERROR RUNNING COMMAND.  Exit code (\d+) from/ ) {
        return $1;
    }

    return 1;
}

1;

