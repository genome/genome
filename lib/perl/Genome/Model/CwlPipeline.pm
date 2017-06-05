package Genome::Model::CwlPipeline;

use strict;
use warnings;

class Genome::Model::CwlPipeline {
    is => 'Genome::Model',
    has_param => {
        main_workflow_file => {
            is => 'Text',
        },
        primary_docker_image => {
            is => 'Text',
            doc => 'docker image for the main toil worker jobs',
        },
    },
};

sub map_workflow_inputs {
    my $self = shift;
    my $build = shift;

    return (build => $build);
}

sub _resolve_workflow_for_build {
    my $self = shift;
    my $build = shift;

    my $wf = Genome::WorkflowBuilder::DAG->create(
        name => $build->workflow_name,
        log_dir => $build->log_directory,
    );

    my $cmd = Genome::WorkflowBuilder::Command->create(
        name => 'Toil Runner',
        command => 'Genome::Model::CwlPipeline::Command::Run',
    );
    $wf->add_operation($cmd);
    $wf->connect_input(
        input_property => 'build',
        destination => $cmd,
        destination_property => 'build'
    );
    $wf->connect_output(
        source => $cmd,
        source_property => 'result',
        output_property => 'result',
    );

    return $wf;
}

1;
