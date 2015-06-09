package Genome::Command::PteroWorkflowMixin;

use strict;
use warnings;

BEGIN {
    # the default 'info' log level is too noisy.
    $ENV{PTERO_PERL_SDK_LOG_LEVEL} = 'warn' unless $ENV{PTERO_PERL_SDK_LOG_LEVEL};
}

use Genome;
use Genome::Utility::Text qw(justify strip_color);
use Try::Tiny qw(try catch);
use Ptero::Proxy::Workflow;
use Ptero::Proxy::Workflow::Execution qw();

my $INDENTATION_STR = '. ';

class Genome::Command::PteroWorkflowMixin {
    is => ['Genome::Command::WithColor'],
    has => [
        summary_threshold => {
            is => 'Int',
            is_optional => 1,
            default_value => 5,
            doc => 'Parallel by stages with this many steps will be summarized. Use zero to summarize every step. Use a negative value to never summarize.'
        },
        logs => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 1,
            doc => 'Display path of log files for running, failed, and errored steps.',
        },
        workflow => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 1,
            doc => 'Display workflow.',
        },
    ],
    has_transient_optional => [
        _executions_of_interest => {
            is => 'ARRAY',
        },
    ],
};

sub _display_ptero_workflow {
    my ($self, $handle, $workflow_name) = @_;

    return unless $self->workflow;
    $handle->say($self->_color_heading("Workflow"));

    my $url = sprintf("%s?name=%s", Genome::Config::get('ptero_workflow_submit_url'),
        $workflow_name);

    my $wf_proxy = try {
        return Ptero::Proxy::Workflow->new($url);
    };
    unless (defined($wf_proxy)) {
        $handle->say("No ptero workflow found at url ($url)");
        return;
    }

    $self->_write_ptero_header($handle);

    $self->_executions_of_interest([]);
    $self->_write_ptero_workflow($handle, $wf_proxy->concrete_workflow, 0, 0);

    if ($self->logs) {
        $self->_write_executions_of_interest($handle);
    }

    return;
}

sub _write_executions_of_interest {
    my $self = shift;
    my $handle = shift;

    $handle->print("\n");
    for my $execution (@{$self->_executions_of_interest}) {
        my $ex_proxy = Ptero::Proxy::Workflow::Execution->new($execution->{details_url});
        if ($ex_proxy->concrete_execution->{status} eq 'errored') {
            $handle->print(join("\n", $self->_color_pair("Name", $ex_proxy->name) . "    " .
                $self->_color_dim("Status: ") . $self->_ptero_status_color($ex_proxy->concrete_execution->{status}),
                $self->_color_pair("Error Message", $ex_proxy->concrete_execution->{data}{errorMessage}),
                $self->_color_pair("Stdout", $ex_proxy->concrete_execution->{data}{stdout}),
                $self->_color_pair("Stderr", $ex_proxy->concrete_execution->{data}{stderr})));
        } else {
            $handle->print(join("\n", $self->_color_pair("Name", $ex_proxy->name) . "    " .
                $self->_color_dim("Status: ") . $self->_ptero_status_color($ex_proxy->concrete_execution->{status}),
                $self->_color_pair("Stderr Log", $ex_proxy->concrete_execution->{data}{stderr_log})));
        }
        $handle->print("\n");
    }
}

my %PTERO_STATUS_COLORS = (
    new => "white",
    scheduled => "white",

    running => "cyan",
    succeeded => "green",

    errored => "red",
    failed => "red",
    canceled => "red",
);

sub _ptero_status_color {
    my ($self, $text) = @_;
    return $self->_colorize_text_by_map($text, $text, %PTERO_STATUS_COLORS);
}

sub _format_line {
    my $self = shift;
    my ($stage, $status, $started, $duration, $pindex, $indent, $name) = @_;

    return join("  ", justify($self->_color_dim($stage), 'right', 8),
        justify($self->_ptero_status_color($status), 'right', 9),
        justify($self->_color_dim("$started"), 'right', 19),
        justify($duration, 'right', 13),
        justify($self->_color_dim($pindex), 'left', 7),
        $INDENTATION_STR x $indent . $name) . "\n";
}

sub _write_ptero_header {
    my $self = shift;
    my $handle = shift;

    $handle->print($self->_color_dim(strip_color($self->_format_line(
        'STAGE',
        'STATUS',
        'STARTED',
        'DURATION',
        'P-INDEX',
        0,
        'NAME'))));

    return;
}

sub _write_ptero_workflow {
    my ($self, $handle, $workflow, $indent, $color) = @_;

    my $execution = $self->_top_level_dag($workflow)->{executions}->{$color};
    if ($execution) {
        $handle->print($self->_format_line(
            '',
            $execution->{status},
            $execution->datetime_started,
            $execution->duration,
            join(', ', $execution->parallel_indexes),
            $indent,
            $workflow->{name}));
    } elsif (scalar(keys %{$workflow->{executions}}) == 0) {
        $handle->print($self->_format_line(
            '',
            $workflow->{status},
            '',
            '',
            '',
            $indent,
            $workflow->{name}));
    }

    my @sorted_tasks = sort {
        $a->{topological_index} <=> $b->{topological_index}}
        $self->_top_level_tasks($workflow);
    for my $task (@sorted_tasks) {
        $self->_write_ptero_task($handle, $task->{name}, $task, $indent+1, $color);
    }

    return;
}

sub _top_level_dag {
    my ($self, $workflow) = @_;

    return $workflow->{tasks}{$workflow->{name}}{methods}[1];
}

sub _top_level_tasks {
    my ($self, $workflow) = @_;

    return values %{$self->_top_level_dag($workflow)->{tasks}};
}

sub _write_ptero_task {
    my $self = shift;
    my ($handle, $task_name, $task, $indent, $color) = @_;

    if ($self->_has_dag_method($task)) {
        $self->_write_ptero_dag_details(@_);
    } else { # has shortcut & execute methods
        $self->_write_ptero_command_details(@_);
    }

    # find executions related to parallel-by split
    my @child_executions = $task->executions_with_parent_color($color);
    return unless scalar(@child_executions);

    if ($self->summary_threshold >= 0 &&
        scalar(@child_executions) >= $self->summary_threshold) {
        $self->_write_ptero_task_summary($handle, $task_name, $task, $indent+1,
            \@child_executions);
    } else {
        for my $child_execution (sort {$a->{color} <=> $b->{color}} @child_executions) {
            $self->_write_ptero_task($handle, $task_name, $task, $indent+1,
                $child_execution->{color});
        }
    }
    return;

}

sub _has_dag_method {
    my $self = shift;
    my $task = shift;

    return (scalar(@{$task->{methods}}) == 1);
}

sub _write_ptero_dag_details {
    my $self = shift;
    my ($handle, $task_name, $task, $indent, $color) = @_;

    my $method = $task->{methods}[0];

    if ($method->{executions}->{$color}) {
        my $execution = $method->{executions}->{$color};
        $handle->print($self->_format_line(
            '',
            $execution->{status},
            $execution->datetime_started,
            $execution->duration,
            join(', ', $execution->parallel_indexes),
            $indent,
            $method->{name}));
    } elsif (scalar(keys %{$method->{executions}}) == 0) {
        $handle->print($self->_format_line(
            '',
            '',
            '',
            '',
            '',
            $indent,
            $method->{name}));
    } else {
        return;
    }

    my @sorted_tasks = sort {
        $a->{topological_index} <=> $b->{topological_index}}
        (values %{$method->{tasks}});
    for my $task (@sorted_tasks) {
        $self->_write_ptero_task($handle, $task->{name}, $task, $indent+1, $color);
    }
}

sub _write_ptero_command_details {
    my $self = shift;
    my ($handle, $task_name, $task, $indent, $color) = @_;

    my ($shortcut, $execute) = @{$task->{methods}};
    if ($self->_method_is_active($shortcut, $color)) {
        $self->_write_ptero_command_details_shortcut(@_);
    } elsif ($self->_method_is_active($execute, $color) ||
             $self->_method_is_failed($execute, $color)) {
        $self->_write_ptero_command_details_execute(@_);
    } else {
        $self->_write_ptero_command_details_unstarted(@_);
    }
}

sub _method_is_active {
    my $self = shift;
    my ($method, $color) = @_;

    my $execution = $method->{executions}->{$color};
    return (defined($execution) && $execution->{status} ne 'failed');
}

sub _method_is_failed {
    my $self = shift;
    my ($method, $color) = @_;

    my $execution = $method->{executions}->{$color};
    return (defined($execution) && $execution->{status} eq 'failed');
}

sub _write_ptero_command_details_shortcut {
    my $self = shift;
    my ($handle, $task_name, $task, $indent, $color) = @_;

    my ($shortcut, $execute) = @{$task->{methods}};
    my $execution = $shortcut->{executions}->{$color};
    if ($execution->{status} eq 'running' ||
        $execution->{status} eq 'failed' ||
        $execution->{status} eq 'errored') {
        push @{$self->_executions_of_interest}, $execution;
    }

    if ($execution) {
        $handle->print($self->_format_line(
            'shortcut',
            $execution->{status},
            $execution->datetime_started,
            $execution->duration,
            join(', ', $execution->parallel_indexes),
            $indent,
            $task_name));
    }
    return;
}

sub _write_ptero_command_details_execute {
    my $self = shift;
    my ($handle, $task_name, $task, $indent, $color) = @_;

    my ($shortcut, $execute) = @{$task->{methods}};
    my $execution = $execute->{executions}->{$color};
    if ($execution->{status} eq 'running' ||
        $execution->{status} eq 'failed' ||
        $execution->{status} eq 'errored') {
        push @{$self->_executions_of_interest}, $execution;
    }

    if ($execution) {
        $handle->print($self->_format_line(
            'execute',
            $execution->{status},
            $execution->datetime_started,
            $execution->duration,
            join(', ', $execution->parallel_indexes),
            $indent,
            $task_name));
    }
    return;
}

sub _write_ptero_command_details_unstarted {
    my $self = shift;
    my ($handle, $task_name, $task, $indent, $color) = @_;

    my $parallel_by_str = '';
    if ($task->{parallel_by}) {
        $parallel_by_str = $self->_color_pair("parallel-by", $task->{parallel_by});
    }

    $handle->printf("%s%s%s\n", justify($parallel_by_str, 'center', 66), $INDENTATION_STR x $indent, $task_name);
    return;
}

sub _write_ptero_task_summary {
    my $self = shift;
    my ($handle, $task_name, $task, $indent, $child_executions) = @_;

    my %statuses;
    for my $execution (@$child_executions) {
        $statuses{$execution->{status}} += 1;
    }
    my $status_str;
    for my $status (sort keys %statuses) {
        $status_str .= sprintf(' %s(%s)', $self->_ptero_status_color($status),
            $statuses{$status});
    }

    $handle->printf("%s%s%s\n", justify($status_str, 'center', 66), $INDENTATION_STR x $indent, $task_name);
    return;
}

1;
