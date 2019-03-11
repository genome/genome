package Genome::Command::CromwellWorkflowMixin;

use strict;
use warnings;

use Genome;
use Genome::Utility::Text qw(justify strip_color);

my $INDENTATION_STR = '. ';
class Genome::Command::CromwellWorkflowMixin {
    is => ['Command::V2'],
    roles => ['Genome::Role::CommandWithColor'],
    has => [
        summary_threshold => {
            is => 'Int',
            is_optional => 1,
            default_value => '5',
            doc => "Scatter stages with this many steps will be summarized. Use zero to summarize every step. Use a negative value to never summarize.",
        },
        logs => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 1,
            doc => 'Display path of log files for running and failed steps.',
        },
        workflow => {
            is => 'Boolean',
            is_optional => 1, 
            default_value => 1,
            doc => 'Display workflow',
        },
    ],
    has_transient_optional => [
        _cromwell_executions_of_interest => {
            is => 'ARRAY',
        },
    ],
};

sub _display_cromwell_workflow {
    my ($self, $handle, $workflow_id) = @_;

    return unless $self->workflow;
    $handle->say($self->_color_heading("Workflow"));

    $self->_cromwell_executions_of_interest([]);

    $self->_write_cromwell_header($handle);
    $self->_write_cromwell_workflow($handle, $workflow_id, 0);

    if ($self->logs) {
        $self->_write_cromwell_executions_of_interest($handle);
    }

    return 1;
}

sub _cromwell_data_for_workflow_id {
    my $self = shift;
    my $workflow_id = shift;

    return  Genome::Cromwell->metadata($workflow_id);
}

sub _write_cromwell_header {
    my $self = shift;
    my $handle = shift;

    $handle->print($self->_color_dim(strip_color($self->_format_cromwell_line(
        qw(ID LSF_ID SHARD STATUS START END NAME), 0
    ))));
}

sub _write_cromwell_workflow {
    my ($self, $handle, $workflow_id, $indent, $shard) = @_;

    my $data = $self->_cromwell_data_for_workflow_id($workflow_id);

    $self->_write_cromwell_workflow_line( $handle, $data, $indent, $shard );

    my @calls = keys (%{ $data->{calls} });
    @calls = sort { $data->{calls}{$a}->[0]{start} cmp $data->{calls}{$b}->[0]{start} } @calls;
    for my $call_name (@calls) {
        for my $call (@{ $data->{calls}{$call_name} }) {
            my $shard = $call->{shardIndex} // -1;

            if ($self->summary_threshold > -1) {
                next if $shard > $self->summary_threshold;

                if ($shard == $self->summary_threshold) {
                    $handle->print($self->_format_cromwell_line('summary threshold reached', '', '', '', '', '', $call_name, $indent+1));
                    next;
                }
            }

            $shard = undef if $shard < 0;

            if (my $subworkflow_id = $call->{subWorkflowId}) {
                $self->_write_cromwell_workflow( $handle, $subworkflow_id, $indent+1, $shard );
            } else {
                $self->_write_cromwell_execution( $handle, $call_name, $call, $indent+1, $shard );
            }
        }
    }
}

sub _write_cromwell_workflow_line {
    my ($self, $handle, $data, $indent, $shard) = @_;

    my $id = $data->{id};
    my $status = $data->{status};
    my $start = $data->{start};
    my $end = $data->{end} // '';
    $shard //= '';
    my $name = $data->{workflowName};
    $handle->print($self->_format_cromwell_line(
        $id, '', $shard, $status, $start, $end, $name, $indent,
    ));

}

sub _write_cromwell_execution {
    my ($self, $handle, $call_name, $call, $indent, $shard) = @_;

    my $job_id = $call->{jobId} // '';
    my $status = $call->{executionStatus};
    my $start = $call->{start};
    my $end = $call->{end} // '';
    $shard //= '';

    $handle->print($self->_format_cromwell_line(
        '', $job_id, $shard, $status, $start, $end, $call_name, $indent,
    ));

    if ($status eq 'Running' or $status eq 'Failed') {
        push @{$self->_cromwell_executions_of_interest},
            [$call_name, $call];
    }

}

sub _format_cromwell_line {
    my $self = shift;
    my ($id, $lsf_job_id, $shard, $status, $start, $end, $name, $indent) = @_;

    return join("  ",
        justify($self->_color_dim($id), 'right', 38),
        justify($lsf_job_id, 'right', 8),
        justify($self->_color_dim($shard), 'left', 7),
        justify($self->_cromwell_status_color($status), 'right', 9),
        justify($self->_color_dim("$start"), 'right', 32),
        justify($end, 'right', 32),
        $INDENTATION_STR x $indent . $name) . "\n";
}

my %CROMWELL_STATUS_COLORS = (
    running => "cyan",
    succeeded => "green",
    done => "green",
    failed => "red",
);

sub _cromwell_status_color {
    my ($self, $text) = @_;
    return $self->_colorize_text_by_map($text, $text, %CROMWELL_STATUS_COLORS);
}

sub _write_cromwell_executions_of_interest {
    my ($self, $handle) = @_;

    return unless $self->logs;

    for my $e (@{ $self->_cromwell_executions_of_interest}) {
        $handle->print(
            $self->_color_pair("Name", $e->[0]) . "    " .
            $self->_color_dim("Status: ") . $self->_cromwell_status_color($e->[1]{executionStatus}),
            "\n",
            $self->_color_pair("Stderr Log", $e->[1]{stderr}),
            "\n",
        );
    }
}

1;
