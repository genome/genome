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
    my %tools = $self->_tools;
    my $process_graph = Genome::WorkflowBuilder::StreamGraph->create(
        output_xml => File::Spec->join($self->temp_staging_directory, "out.xml"),
    );
    my %process_ref;
    while (my ($name, $tool) = each %tools) {
        $process_ref{$name} = Genome::WorkflowBuilder::StreamProcess->create(
            name => $name,
            args => [$tool->cmd_line],
            in_file_link => $self->_input_file_for_tool($tool, $name),
            out_file_link => $self->_output_file_for_tool($name),
            err_file_link => $self->_error_file_for_tool($name),
        );
        $process_graph->add_process($process_ref{$name});
    }

    while (my ($name, $tool) = each %tools) {
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

    for my $tool (values %tools) {
        $self->_add_metrics($tool);
    }
    return 1;
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
    my ($name, $gmt_params) = @_;
    if (defined $name->output_file_accessor) {
        my $output_param_name = $name->output_file_accessor;
        $gmt_params->{$output_param_name} = Genome::Sys->create_temp_file;
    }
    return $name->create(gmt_params => $gmt_params);
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

