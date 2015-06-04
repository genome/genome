package Genome::Ptero::Wrapper;

use strict;
use warnings FATAL => qw(all);
use Genome;
use Data::Dump qw();
use IO::Handle;
use Cwd qw(abs_path);
use Genome::Utility::Text qw(
    sanitize_string_for_filesystem
);
use Genome::Utility::Inputs qw(encode decode);
use Data::Dump qw(pp);
use Ptero::Proxy::Workflow::Execution;
use Try::Tiny qw(try catch);

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
        'log_directory' => {
            is => 'Path',
            doc => "Where the log files should be stored (Their names will " .
                "be automatically generated.)",
        },
    ],
    has_transient_optional => [
        'execution' => {
            is => 'Ptero::Proxy::Workflow::Execution',
        }
    ],

    doc => 'Ptero execution wrapper for Genome Command objects',
};


sub execute {
    my $self = shift;

    $self->debug_message("Creating log directory %s", $self->log_directory);
    Genome::Sys->create_directory($self->log_directory);
    validate_environment();

    $self->execution(Ptero::Proxy::Workflow::Execution->new(
        $ENV{PTERO_WORKFLOW_EXECUTION_URL}));

    # this will get saved to the execution's stderr/stdout
    $self->_log_execution_information;

    $self->_setup_logging;

    # this will get logged to the log files
    $self->_log_execution_information;

    my ($status_ok, $command) = $self->_instantiate_and_run_command;

    $self->_teardown_logging;

    unless ($status_ok) {
        Carp::croak sprintf("Failed to instantiate and run (%) on (%s)\n",
            $self->method, $self->command_class);
    }

    printf SAVED_STDERR "Setting outputs: %s\n", pp(_get_command_outputs($command, $self->command_class));
    $self->execution->set_outputs(
        _get_command_outputs($command, $self->command_class));

    return 1;
}

sub validate_environment {
    my $self = shift;
    unless (defined $ENV{PTERO_WORKFLOW_EXECUTION_URL}) {
        die $self->error_message("Environment variable PTERO_WORKFLOW_EXECUTION_URL must be set");
        exit 1;
    }
}

sub _stdout_log_path {
    my $self = shift;

    my $base_name = $self->execution->name;
    my $output_log = File::Spec->join(abs_path($self->log_directory),
        sanitize_string_for_filesystem("$base_name.out"));
}

sub _stderr_log_path {
    my $self = shift;

    my $base_name = $self->execution->name;
    my $output_log = File::Spec->join(abs_path($self->log_directory),
        sanitize_string_for_filesystem("$base_name.err"));
}

sub _setup_logging {
    my $self = shift;

    $self->debug_message(
        "Preparing to redirect stderr to (%s) and stdout to (%s)",
        $self->_stderr_log_path, $self->_stdout_log_path
    );

    $self->execution->update_data(
        stdout_log => $self->_stdout_log_path,
        stderr_log => $self->_stderr_log_path,
    );

    open(SAVED_STDOUT, ">&STDOUT") || die "Can't save STDOUT\n";
    open(SAVED_STDERR, ">&STDERR") || die "Can't save STDERR\n";

    # redirect stdout/stderr through the annotate-log filter
    open OUTPUT, '|-',  'annotate-log cat > ' .
        $self->_stdout_log_path or die $!;
    open ERROR, '|-',  'annotate-log cat > ' .
        $self->_stderr_log_path or die $!;

    STDOUT->fdopen(\*OUTPUT, 'w') or die $!;
    STDERR->fdopen(\*ERROR,  'w') or die $!;
}

sub _teardown_logging {
    my $self = shift;

    printf SAVED_STDERR "Removing stderr and stdout redirection\n";

    open(STDOUT, ">&SAVED_STDOUT") || die "Can't restore STDOUT\n";
    open(STDERR, ">&SAVED_STDERR") || die "Can't restore STDERR\n";
}

sub _log_execution_information {
    my $self = shift;

    $self->status_message("COMMAND: PTERO_WORKFLOW_EXECUTION_URL=%s %s " .
        "ptero wrapper --command-class=\"%s\" --method=\"%s\" --log-directory=\"%s\"",
        $ENV{PTERO_WORKFLOW_EXECUTION_URL}, $0,
        $self->command_class, $self->method, abs_path($self->log_directory));
}

sub _instantiate_and_run_command {
    my $self = shift;

    printf SAVED_STDERR "Instantiating command %s\n", $self->command_class;

    my $cmd;
    my $status_ok = $self->_eval_command_class;
    $status_ok = defined( $cmd = $self->_instantiate_command )  if $status_ok;
    $status_ok = $self->_run_command($cmd)                      if $status_ok;
    $status_ok = $self->_commit                                 if $status_ok;

    if ($status_ok) {
        $self->status_message("Succeeded to %s command %s",
            $self->method, $self->command_class)
    }

    return ($status_ok, $cmd);
}

sub _eval_command_class {
    my $self = shift;

    my $eval_succeeded = try {
        my $pkg = $self->command_class;
        eval "use $pkg";
        return 1;
    }
    catch {
        Carp::cluck sprintf("Instantiating class (%s) failed\n", $self->command_class);
        return;
    };

    return $eval_succeeded;
}

sub _instantiate_command {
    my $self = shift;
    my $inputs;

    my $pkg = $self->command_class;
    my $cmd = try {
        $inputs = decode($self->execution->inputs);
        return $pkg->create(%$inputs);
    } catch {
        if (defined $inputs) {
            Carp::cluck sprintf(
                "Failed to instantiate class (%s) with inputs (%s): %s",
                $pkg, Data::Dump::pp($inputs), $_)
        }
        else {
            Carp::cluck sprintf(
                "Failed to instantiate class (%s) with : %s",
                $pkg, Data::Dump::pp($inputs), $_)
        }
        return;
    };

    return $cmd;
}

sub _run_command {
    my ($self, $command) = @_;

    printf SAVED_STDERR "Running command %s\n", $self->command_class;

    my $method = $self->method;
    my $ret = try {
        return $command->$method();
    } catch {
        Carp::cluck sprintf(
            "Crashed in %s for command %s: %s",
            $self->method, $self->command_class, $_,
        );
        return;
    };

    unless ($ret) {
        Carp::cluck sprintf("Failed to %s for command %s.",
            $self->method, $self->command_class,
        );
    }

    return $ret;
}

sub _commit {
    my $self = shift;

    my $rv = try {
        return UR::Context->commit();
    } catch {
        Carp::confess "Failed to commit: $_";
        return;
    };

    unless ($rv) {
        Carp::confess "Failed to commit: see previously logged errors";
    }
    return $rv;
}

sub _get_command_outputs {
    my ($cmd, $pkg) = @_;

    my %outputs;
    for my $prop (_output_properties($pkg)) {
        my $prop_name = $prop->property_name;
        my $value = $prop->is_many ? [$cmd->$prop_name] : $cmd->$prop_name;
        $outputs{$prop_name} = $value;
    }

    return encode(\%outputs);
}

sub _output_properties {
    my $pkg = shift;

    return $pkg->__meta__->properties(is_output => 1);
}


1;
