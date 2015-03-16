package Genome::Qc::Result;

use strict;
use warnings;
use Genome;
use Genome::Qc::Factory;
use IPC::Run qw(run);

class Genome::Qc::Result {
    is => 'Genome::SoftwareResult::StageableSimple',
    has_input => [
        alignment_result => {
            is => 'Genome::InstrumentData::AlignedBamResult',
        },
        #This will need to be backfilled once we have a db-backed Qc::Config
        #qc_config => {
        #    is => 'Genome::Qc::Config',
        #},
    ],
};

sub qc_config {
    return Genome::Qc::Factory::get_config(name => "default");
}

sub get_metrics {
    my $self = shift;
    my @metrics = Genome::SoftwareResult::Metric->get(
        software_result => $self,
    );
    return map {$_->metric_name => $_->metric_value} @metrics;
}

sub _run {
    my $self = shift;
    $self->_run_tools;
    my $temp = File::Spec->join($self->temp_staging_directory, "temp");
    `touch $temp`;
    return 1;
}

sub _run_tools {
    my $self = shift;

    if ($self->_streaming_tools) {
        $self->_run_streaming_tools;
    }

    my %tools = $self->_tools;
    for my $tool (map {$tools{$_}} $self->_non_streaming_tools) {
        $self->_run_tool_simple($tool);
    }

    for my $tool (map {$tools{$_}} $self->_streaming_tools, $self->_non_streaming_tools) {
        $self->_add_metrics($tool);
    }
}

sub _run_streaming_tools {
    my $self = shift;
    my %tools = $self->_tools;
    my $process_graph = Genome::WorkflowBuilder::StreamGraph->create();
    my %process_ref;
    for my $name ($self->_streaming_tools) {
        my $tool = $tools{$name};
        $process_ref{$name} = Genome::WorkflowBuilder::StreamProcess->create(
            name => $name,
            args => [$tool->cmd_line],
            in_file_link => $self->_input_file_for_tool($tool, $name),
            out_file_link => $self->_output_file_for_tool($name),
            err_file_link => $self->_error_file_for_tool($name)
        );
        $process_graph->add_process($process_ref{$name});
    }

    for my $name ($self->_streaming_tools) {
        my $tool = $tools{$name};
        my $dependency = $self->_dependency_for_tool($name);
        if (defined $dependency) {
            my $link = Genome::WorkflowBuilder::StreamLink->create(
                source => $process_ref{$dependency->{name}},
                target => $process_ref{$name},
                source_fd => $dependency->{fd},
            );
            $process_graph->add_link($link);
        }
    }
    $process_graph->execute;
}

sub _dependency_for_tool {
    my ($self, $name) = @_;
    return $self->qc_config->get_commands_for_alignment_result->{$name}->{dependency};
}

sub _input_file_for_tool {
    my ($self, $tool, $name) = @_;
    my $input_file_method = $self->qc_config->get_commands_for_alignment_result->{$name}->{in_file};
    if (defined $input_file_method) {
        return $tool->$input_file_method;
    }
    return undef;
}

sub _error_file_for_tool {
    my ($self, $name) = @_;
    my $file_name = $self->qc_config->get_commands_for_alignment_result->{$name}->{error_file};
    if (defined $file_name) {
        return File::Spec->join($self->temp_staging_directory, $file_name);
    }
    return undef;
}

sub _output_file_for_tool {
    my ($self, $name) = @_;
    my $file_name = $self->qc_config->get_commands_for_alignment_result->{$name}->{out_file};
    if (defined $file_name) {
        return File::Spec->join($self->temp_staging_directory, $file_name);
    }
    return undef;
}

sub _run_tool_simple {
    my $self = shift;
    my $tool = shift;
    run($tool->cmd_line);
}

sub _tools {
    my $self = shift;
    my $commands = $self->qc_config->get_commands_for_alignment_result($self->alignment_result);
    my %tools;
    for my $name (keys %$commands) {
        my $tool = _tool_from_name_and_params($commands->{$name}->{class},
                    $commands->{$name}->{params});
        $tools{$name} = $tool;
    }
    return %tools;
}

sub _streaming_tools {
    my $self = shift;
    my %tools = $self->_tools;
    return grep {$tools{$_}->supports_streaming} keys %tools;
}

sub _non_streaming_tools {
    my $self = shift;
    my %tools = $self->_tools;
    return grep {!$tools{$_}->supports_streaming} keys %tools;
}

sub _tool_from_name_and_params {
    my ($name, $params) = @_;

    return $name->create(%{$params});
}

sub _add_metrics {
    my ($self, $tool) = @_;
    my %metrics = %{$tool->get_metrics};
    while (my ($name, $value) = each %metrics) {
        $self->add_metric(
            metric_name => $name,
            metric_value => $value,
        );
    }
}

1;

