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
use Try::Tiny qw(try catch finally);

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

    my $command = try {
        my $cmd = $self->_instantiate_command(decode($self->execution->inputs));
        $self->_run_command($cmd);
        return $cmd;
    } catch {
        die $_;
    } finally {
        $self->_teardown_logging;
    };

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

sub _instantiate_command {
    my ($self, $inputs) = @_;

    printf SAVED_STDERR "Instantiating command %s\n", $self->command_class;

    my $pkg = $self->command_class;
    my $cmd = try {
        eval "use $pkg";
        $pkg->create(%$inputs)
    } catch {
        Carp::confess sprintf(
            "Failed to instantiate class (%s) with inputs (%s): %s",
            $pkg, Data::Dump::pp($inputs), $_)
    };

    return $cmd;
}

sub _run_command {
    my ($self, $command) = @_;

    printf SAVED_STDERR "Running command %s\n", $self->command_class;

    my $method = $self->method;
    my $ret = try {
        $command->$method()
    } catch {
        Carp::confess sprintf(
            "Crashed in %s for command %s: %s",
            $self->method, $self->command_class, $_,
        );
    };
    unless ($ret) {
        Carp::confess sprintf("Failed to %s for command %s.",
            $self->method, $self->command_class,
        );
    }

    $self->_commit();

    $self->status_message("Succeeded to %s command %s",
        $method, $self->command_class);
}

sub _commit {
    my $self = shift;
    my $rv = try {
        UR::Context->commit()
    } catch {
        Carp::confess "Failed to commit: $_";
    };

    unless ($rv) {
        Carp::confess "Failed to commit: see previously logged errors";
    }
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
