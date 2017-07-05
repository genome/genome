package Genome::Model::CwlPipeline;

use strict;
use warnings;

use feature qw(switch);
use Genome;

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


#handle specification of inputs such as might be found in an AnP Config YAML
sub process_input_data {
    my $self = shift;
    my %data = @_;

    for my $name (keys %data) {
        my $value = $self->determine_input_object($name, $data{$name});
        unless ($value) {
            $self->fatal_message('Failed to determine input for name "%s" and value "%s".', $name, $data{$name});
        }

        $self->add_input(
            name => $name,
            value_id => $value->id,
            value_class_name => $value->class,
        );
    }

    return 1;
}

sub determine_input_object {
    my $self = shift;
    my $name = shift;
    my $value_identifier = shift;

    #This could use something like Command::Dispatch::Shell->resolve_param_value_from_text
    #to allow for non-ID values
    for ($name) {
        when ('instrument_data')          {
            return Genome::InstrumentData->get($value_identifier);
        }
        when ('reference_build')          {
            return Genome::Model::Build::ReferenceSequence->get($value_identifier);
        }
        when ('annotation_build')         {
            return Genome::Model::Build::ImportedAnnotation->get($value_identifier);
        }
        default                           {
            return UR::Value::Text->get($value_identifier);
        }
    }

    return;
}


1;
