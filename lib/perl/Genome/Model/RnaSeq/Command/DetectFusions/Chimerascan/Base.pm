package Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::Base;

use strict;
use warnings;

use Genome;

class Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::Base {
    is => 'Genome::Model::RnaSeq::Command::DetectFusions::Base',
};

sub command_class_prefix {
    die "Abstract";
}

sub execute {
    my $self = shift;

    my $dag = $self->build_workflow();
    my %inputs = $self->build->model->map_workflow_inputs;
    my $outputs = $dag->execute(%inputs);

    $self->result($outputs->{'result'});
    return 1;
}

sub build_workflow {
    my $self = shift;

    my $dag = Genome::WorkflowBuilder::DAG->create(
        name => 'ChimerascanWorkflow',
        log_dir => $self->build->log_directory,
    );

    my $index_command = Genome::WorkflowBuilder::Command->create(
        name => 'IndexCommand',
        command => $self->index_command_name,
    );
    $dag->add_operation($index_command);

    my $chimerascan_command = Genome::WorkflowBuilder::Command->create(
        name => 'ChimerascanCommand',
        command => $self->chimerascan_command_name,
    );
    $dag->add_operation($chimerascan_command);

    $self->_add_inputs($dag, $index_command);
    $self->_add_inputs($dag, $chimerascan_command);

    $self->_add_links($dag, $index_command, $chimerascan_command);
    $self->_add_outputs($dag, $chimerascan_command);

    return $dag;
}

sub _add_inputs {
    my $self = shift;
    my $dag = shift;
    my $command = shift;

    $dag->connect_input(
        input_property => 'detector_version',
        destination => $command,
        destination_property => 'detector_version',
    );

    $dag->connect_input(
        input_property => 'detector_params',
        destination => $command,
        destination_property => 'detector_params',
    );

    $dag->connect_input(
        input_property => 'build',
        destination => $command,
        destination_property => 'build',
    );
}

sub _add_links {
    my $self = shift;
    my $source_op = shift;
    my $destination_op = shift;

    my $link = Genome::WorkflowBuilder::Link->create(
        source => $source_op, source_property => 'single_output',
        destination => $destination_op, destination_property => 'input',
    );
}

sub _add_outputs {
    my $self = shift;
    my $dag = shift;
    my $command = shift;

    $dag->connect_input(
        input_property => 'build_id',
        destination => $command,
        destination_property => 'build_id',
    );
}

sub chimerascan_command_name {
    my $self = shift;

    return $self->command_class_prefix . "::Detector";
}

sub index_command_name {
    my $self = shift;

    return $self->command_class_prefix . "::Index";
}


