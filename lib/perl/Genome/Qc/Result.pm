package Genome::Qc::Result;

use strict;
use warnings;
use Genome;
use Genome::Qc::Factory;

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
    my $process_graph = Genome::Sys::ProcessGraph->create;
    my %tools = $self->_streaming_tools;

    while (my ($name, $tool)  = each %tools) {
        $process_ref{$name} = $process_graph->add(
            command_line => $tool->cmd_line,
            stderr => ?);
    }

    while (my ($name, $tool) = each %tools) {
        my $dependency = $self->_dependency_for_tool($name);
        if (defined $dependency) {
            $process_graph->connect($process_ref{$dependency->{name}},
                $dependency->{fd}, $process_ref{$name}, "STDIN");
        }
    }

    $process_graph->execute;

    for my $tool ($self->_non_streaming_tools) {
        $self->_run_tools_simple($tool);
    }
}

sub _dependency_for_tool {
    my ($self, $name) = @_;
    return $self->qc_config->get_commands_for_alignment_results->{$name}->{dependency};
}

sub _run_tools_simple {
    my $self = shift;
    my $tool = shift;
    $tool->run;
    $self->_add_metrics($tool);
}

sub _tools {
    my $self = shift;
    my $commands = $self->qc_config->get_commands_for_alignment_result($self->alignment_result);
    my @tools;
    for my $name (keys %$commands) {
        my $tool = _tool_from_name_and_params($name, $commands->{$name}->{params});
        push @tools, $tool;
    }
    return @tools;
}

sub _streaming_tools {
    return grep {$_->supports_streaming} $self->_tools;
}

sub _non_streaming_tools {
    return grep {!$_->supports_streaming} $self->_tools;
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

