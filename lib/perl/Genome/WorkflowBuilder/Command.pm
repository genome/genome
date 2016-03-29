package Genome::WorkflowBuilder::Command;

use strict;
use warnings;

use Genome;
use Cwd qw();
use Genome::Sys::LSF::ResourceParser qw(parse_lsf_params);
use Data::UUID;
use Carp qw();
use Data::Dump qw(pp);
use Try::Tiny;


class Genome::WorkflowBuilder::Command {
    is => 'Genome::WorkflowBuilder::Detail::Operation',

    has => [
        command => {
            is => 'Command',
        },
    ],
    has_optional => [
        lsf_queue => {
            is => 'String',
        },
        lsf_project => {
            is => 'String',
        },
        lsf_resource => {
            is => 'String',
        },
    ],
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    eval sprintf("require %s", $self->command);
    my $error = $@;
    if ($error) {
        Carp::confess(sprintf("Failed to load command class (%s): %s\n",
                $self->command, $error));
    }
    return $self;
}

sub get_ptero_builder_task {
    require Ptero::Builder::Detail::Workflow::Task;

    my $self = shift;
    my $log_dir = shift;

    $self->validate;

    my %params = (
        name => $self->name,
        methods => [
            $self->_get_ptero_shortcut_method($log_dir),
            $self->_get_ptero_execute_method($log_dir),
        ],
    );
    if (defined $self->parallel_by) {
        $params{parallel_by} = $self->parallel_by;
    }
    return Ptero::Builder::Detail::Workflow::Task->new(%params);
}

sub _get_ptero_shortcut_method {
    require Ptero::Builder::Job;

    my $self = shift;
    my $log_dir = shift;

    my $retry_exit_code = Genome::Config::get('software_result_async_locking_exit_code');
    my $initial_interval = Genome::Config::get('software_result_async_locking_initial_interval');
    my $max_interval = Genome::Config::get('software_result_async_locking_max_interval');
    my $num_attempts = Genome::Config::get('software_result_async_locking_num_attempts');
    my $retry_settings = {
        exitCode => $retry_exit_code + 0,
        initialInterval => $initial_interval + 0,
        maxInterval => $max_interval + 0,
        attempts => $num_attempts + 0,
    };

    # Set env variable to allow genome to exit instead of wait while
    # attempting shortcutting
    my $scope_guard = Genome::Config::set_env("software_result_async_locking", 1);
    my $env = $self->_get_sanitized_env();

    my %job_args = (
        name => 'shortcut',
        service_url => Genome::Config::get('ptero_shell_command_service_url'),
        service_data_to_save => ['error_message', 'error'],
        parameters => {
            commandLine => [
                'genome', 'ptero', 'wrapper',
                '--command-class', $self->command,
                '--method', 'shortcut',
                '--log-directory', $log_dir,
            ],
            environment => $env,
            retrySettings => $retry_settings,
            user => Genome::Sys->username,
            workingDirectory => Cwd::getcwd,
        },
    );


    return Ptero::Builder::Job->new(%job_args);
}

sub _get_ptero_execute_method {
    require Ptero::Builder::Job;


    my $self = shift;
    my $log_dir = shift;
    my $ptero_lsf_parameters = $self->_get_ptero_lsf_parameters();
    $ptero_lsf_parameters->{command} = sprintf(
        'genome ptero wrapper --command-class %s '
        .'--method execute --log-directory %s',
        $self->command, $log_dir);
    $ptero_lsf_parameters->{environment} = $self->_get_sanitized_env();
    $ptero_lsf_parameters->{user} = Genome::Sys->username;
    $ptero_lsf_parameters->{cwd} = Cwd::getcwd;
    $ptero_lsf_parameters->{pollingInterval} =
        Genome::Config::get('ptero_lsf_polling_interval') + 0;

    my $project_name = Genome::Config::get('lsf_project_name');
    if ($project_name) {
        $ptero_lsf_parameters->{options}{projectName} = $project_name;
    }

    return Ptero::Builder::Job->new(
        name => 'execute',
        service_url => Genome::Config::get('ptero_lsf_service_url'),
        service_data_to_save => ['error_message', 'error', 'lsfJobId'],
        parameters => $ptero_lsf_parameters);
}

sub _get_ptero_lsf_parameters {
    my $self = shift;
    my %attributes = $self->operation_type_attributes;
    my $lsf_params = parse_lsf_params( $attributes{lsfResource} );

    my $set_lsf_option = sub {
        my ($option, $value) = @_;
        if (defined($value) and length($value) and
            not exists($lsf_params->{options}->{$option})) {
            $lsf_params->{options}->{$option} = $value;
        }
    };

    $set_lsf_option->('queue', $attributes{lsfQueue});
    $set_lsf_option->('projectName', $attributes{lsfProject});

    my $default_job_group = join('/',
        Genome::Config::get('lsf_job_group'),
        Genome::Sys->username,
    );
    $set_lsf_option->('jobGroup', $default_job_group);

    my ($stderr, $stdout, $postexec) = _get_lsf_log_paths();
    $lsf_params->{options}->{errFile} = $stderr;
    $lsf_params->{options}->{outFile} = $stdout;

    # Redirect stdout/err of post-exec command to file for debugging.
    # Clean up that file if post-exec exits normally.
    $lsf_params->{options}->{postExecCmd} = sprintf(
        "bash -c '%s' > %s 2>&1 && rm -f %s",
        "ptero-lsf-post-exec --stderr $stderr --stdout $stdout",
        $postexec, $postexec
    );

    # Unlike post-exec, stdout/err of pre-exec command gets written to
    # errFile or outFile by lsf.  We want to exit 0 no matter what so
    # that the job doesn't get caught in a PEND->RUN->PEND loop.
    $lsf_params->{options}->{preExecCmd} = "ptero-lsf-pre-exec; exit 0;";

    return $lsf_params;
}

sub _get_lsf_log_paths {
    my $ug = Data::UUID->new();
    my $uuid = $ug->create();
    my $uuid_str = $ug->to_string($uuid);

    my $base = sprintf("/tmp/ptero-lsf-logfile-%s", $uuid_str);
    return ($base . '.err', $base . '.out', $base . "-postexec.log");
}


# ------------------------------------------------------------------------------
# Inherited Methods
# ------------------------------------------------------------------------------

sub from_xml_element {
    my ($class, $element) = @_;

    return $class->create(
        name => $element->getAttribute('name'),
        parallel_by => $element->getAttribute('parallelBy'),
        $class->operationtype_attributes_from_xml_element($element),
    );
}

sub expected_attributes {
    return (
        command => 'commandClass',
        lsf_project => 'lsfProject',
        lsf_queue => 'lsfQueue',
        lsf_resource => 'lsfResource',
    );
}

sub input_properties {
    my $self = shift;

    my @metas = $self->command->__meta__->properties(
        is_input => 1, is_optional => 0);

    my @metas_without_defaults = grep {! defined($_->default_value)} @metas;

    my @result = map {$_->property_name} @metas_without_defaults;
    return sort @result;
}

sub operation_type_attributes {
    my $self = shift;
    my %attributes = (
        commandClass => $self->command,
    );
    my %expected_attributes = $self->expected_attributes;
    for my $name (keys(%expected_attributes)) {
        my $value;
        if (defined($self->$name)) {
            $value = $self->$name;
        } else {
            $value = $self->_get_attribute_from_command($name);
        }

        if (defined($value)) {
            $attributes{$expected_attributes{$name}} = $value;
        }
    }
    return %attributes;
}

sub output_properties {
    my $self = shift;
    return sort map {$_->property_name} $self->command->__meta__->properties(
        is_output => 1);
}

sub validate {
    my $self = shift;

    $self->SUPER::validate();

    if (defined($self->parallel_by)) {
        my $prop = $self->command->__meta__->properties(
            property_name => $self->parallel_by, is_input => 1);
        if (!defined($prop)) {
            die sprintf("Failed to verify that requested " .
                    "parallel_by property '%s' was an input " .
                    "on command (%s)",
                    $self->parallel_by, $self->command->class);
        }
    }
}

sub is_input_property {
    my ($self, $property_name) = @_;
    return $self->command->__meta__->properties(
            property_name => $property_name, is_input => 1)
        || $self->command->__meta__->properties(
            property_name => $property_name, is_param => 1);
}

sub is_output_property {
    my ($self, $property_name) = @_;
    return $self->command->__meta__->properties(property_name => $property_name,
        is_output => 1);
}

sub is_many_property {
    my ($self, $property_name) = @_;
    return $self->command->__meta__->properties(property_name => $property_name,
        is_many => 1);
}


# ------------------------------------------------------------------------------
# Private Methods
# ------------------------------------------------------------------------------

sub _get_attribute_from_command {
    my ($self, $property_name) = @_;

    my $property = $self->command->__meta__->properties(
        property_name => $property_name);
    return unless defined $property;

    if (defined $property->default_value) {
        return $property->default_value;
    }
    elsif ($property->calculated_default) {
        return $property->calculated_default->();
    }
    else {
        return;
    }
}

sub _execute_inline {
    my ($self, $inputs) = @_;

    my $cmd = $self->_instantiate_command($inputs);
    $self->_run_command($cmd);
    return _get_command_outputs($cmd, $self->command);
}

sub _instantiate_command {
    my ($self, $inputs) = @_;

    $self->status_message("Instantiating command %s", $self->command);

    my $pkg = $self->command;
    my $cmd = try {
        eval "use $pkg";
        $pkg->create(%$inputs)
    } catch {
        Carp::confess sprintf(
            "Failed to instantiate class (%s) with inputs (%s): %s",
            $pkg, pp($inputs), pp($_))
    };

    return $cmd;
}

sub _run_command {
    my ($self, $cmd) = @_;

    $self->status_message("Running command %s", $self->command);

    my $ret = try {
        $cmd->execute()
    } catch {
        Carp::confess sprintf(
            "Crashed in execute for command %s: %s",
            $self->command, $_,
        );
    };
    unless ($ret) {
        Carp::confess sprintf("Failed to execute for command %s.",
            $self->command,
        );
    }

    $self->status_message("Succeeded to execute command %s", $self->command);
}

sub _get_command_outputs {
    my ($cmd, $pkg) = @_;

    my %outputs;
    my @output_properties = $pkg->__meta__->properties(is_output => 1);
    for my $prop (@output_properties) {
        my $prop_name = $prop->property_name;
        my $value = $prop->is_many ? [$cmd->$prop_name] : $cmd->$prop_name;
        $outputs{$prop_name} = $value;
    }

    return \%outputs;
}


1;
