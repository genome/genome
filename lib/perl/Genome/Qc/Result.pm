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
    my %tools = $self->_tools;
=cut
    my $process_graph = Genome::Sys::ProcessGraph->create(
        log_directory => File::Spec->join($self->temp_staging_directory, "process_logs"));
    my %process_ref;
    for my $name ($self->_streaming_tools) {
        my $tool = $tools{$name};
        my $log_file = File::Spec->join($self->temp_staging_directory, "$name.log");
        $process_ref{$name} = $process_graph->add(
            command_line => $tool->cmd_line,
            stderr => $log_file);
    }

    for my $name ($self->_streaming_tools) {
        my $tool = $tools{$name};
        my $dependency = $self->_dependency_for_tool($name);
        if (defined $dependency) {
            $process_graph->connect($process_ref{$dependency->{name}},
                $dependency->{fd}, $process_ref{$name}, "STDIN");
        }
        my $in_file = $self->_input_file_for_tool($name);
        if (defined $in_file) {
            $process_graph->connect_input_file($tool->$in_file, $process_ref{$name}, "STDIN");
        }
        my $out_file = File::Spec->join($self->temp_staging_directory, $self->_output_file_for_tool($name));
        if (defined $out_file) {
            $process_graph->connect_output_file($tool->$out_file, $process_ref{$name}, "STDOUT");
        }
        my $err_file = File::Spec->join($self->temp_staging_directory, $self->_error_file_for_tool($name));
        if (defined $err_file) {
            $process_graph->connect_output_file($tool->$err_file, $process_ref{$name}, "STDERR");
        }
    }

    $process_graph->execute;
=cut
    for my $tool (map {$tools{$_}} $self->_non_streaming_tools) {
        $self->_run_tool_simple($tool);
    }

    for my $tool (map {$tools{$_}} $self->_streaming_tools, $self->_non_streaming_tools) {
        $self->_add_metrics($tool);
    }
}

sub _dependency_for_tool {
    my ($self, $name) = @_;
    return $self->qc_config->get_commands_for_alignment_results->{$name}->{dependency};
}

sub _input_file_for_tool {
    my ($self, $name) = @_;
    return $self->qc_config->get_commands_for_alignment_results->{$name}->{in_file};
}

sub _error_file_for_tool {
    my ($self, $name) = @_;
    return $self->qc_config->get_commands_for_alignment_results->{$name}->{error_file};
}

sub _output_file_for_tool {
    my ($self, $name) = @_;
    return $self->qc_config->get_commands_for_alignment_results->{$name}->{out_file};
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

