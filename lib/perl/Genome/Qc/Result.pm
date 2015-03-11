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
    for my $tool ($self->_tools) {
        $tool->run;
        $self->_add_metrics($tool);
    }
    my $temp = File::Spec->join($self->temp_staging_directory, "temp");
    `touch $temp`;
    return 1;
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

sub _tool_from_name_and_params {
    my ($name, $params) = @_;

    my $tool = "Genome::Qc::Tool::$name";
    return $tool->create(%{$params});
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

