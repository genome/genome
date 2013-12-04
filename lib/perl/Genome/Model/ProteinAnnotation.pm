package Genome::Model::ProteinAnnotation;

use strict;
use warnings;
use Genome;

class Genome::Model::ProteinAnnotation {
    is => 'Genome::Model',
    has => [
        subject => {
            is => 'Genome::Taxon',
            id_by => 'subject_id',
            doc => 'taxon of species from which input gene prediction originated',
        },
    ],
    has_param => [
        strategy => {
            is => 'Text',
            doc => 'a string describing what predictors should be run and what should be done with their results',
        },
        chunk_size => {
            is => 'Number',
            doc => 'number of max sequences per fasta chunk should the fasta need to be divided up',
        },
        dump_predictions_to_file => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 1,
        },
        dump_predictions_to_biosql => {
            is => 'Boolean',
            is_optional => 1,
            default => 1,
        },
    ],
    has_input => [
        input_fasta_file => {
            is => 'Text',
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'input_fasta_file', value_class_name => 'UR::Value::Text' ], 
            is_mutable => 1,
        },
    ],
    has_transient => [
        _workflow_inputs => {
            is_optional => 1,
        },
    ],
    doc => 'execute protein prediction and uploads results to a database',
};

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__;

    unless(-e $self->input_fasta_file && -f $self->input_fasta_file) {
        push @errors,UR::Object::Tag->create(
            type => 'error',
            properties => ['input_fasta_file'],
            desc => 'input_fasta_file does not exist or is not a file'
        );
    }

    return @errors;
}

sub server_dispatch {
    return $ENV{GENOME_LSF_QUEUE_BUILD_WORKER};
}

sub map_workflow_inputs {
    my ($self, $build) = @_;
    unless ($self->_workflow_inputs) {
        my %workflow_inputs;
        for my $type (qw/ is_input is_param /) {
            my @properties = $self->__meta__->properties($type => 1);
            for my $property (@properties) {
                my $name = $property->property_name;
                $workflow_inputs{$name} = $self->$name;
            }
        }
        
        # Manually add a few things
        $workflow_inputs{output_directory} = $build->data_directory;
        $workflow_inputs{gram_stain} = $self->subject->gram_stain_category;
        $workflow_inputs{dump_workflow_png_file} = 1;
        $workflow_inputs{dump_workflow_xml_file} = 1;
        $self->_workflow_inputs(\%workflow_inputs);
    }
    
    my %inputs = %{$self->_workflow_inputs};
    return map { $_ => $inputs{$_} } sort keys %inputs;
}

sub _resolve_workflow_for_build {
    my ($self, $build, $lsf_queue, $lsf_project) = @_;
    my %inputs = $self->map_workflow_inputs($build);
    if (defined $lsf_queue and $lsf_queue eq 'inline') {
        $inputs{run_inline} = 1;
    }
    else {
        $inputs{run_inline} = 0;
    }

    my $workflow = Workflow::Model->create(
        name => $build->workflow_name,
        input_properties => [ sort keys %inputs ],
        output_properties => ['bio_seq_features'],
        log_dir => $build->log_directory,
    );

    # Generate prediction workflow
    my $prediction_operation = $self->_generate_prediction_operation($workflow, %inputs);
    my $db_upload_operation = $self->_generate_upload_operation($workflow, $prediction_operation, %inputs);
    return $workflow;
}

sub _generate_prediction_operation {
    my ($self, $workflow, %inputs) = @_;

    my $prediction_operation = $workflow->add_operation(
        operation_type => Workflow::OperationType::Command->create(
            command_class_name => 'Genome::Model::Tools::Predictor::WorkflowGenerator',
        ),
        name => 'prediction workflow',
    );

    my @prediction_input_properties = Genome::Model::Tools::Predictor::WorkflowGenerator->__meta__->properties(is_input => 1);
    for my $property (@prediction_input_properties) {
        my $name = $property->property_name;
        unless (exists $inputs{$name} and defined $inputs{$name}) {
            if ($property->is_optional) {
                next; # Property is optional, so we don't need to worry too much about linking it
            }
            else {
                die "Prediction workflow requires an input $name, but " . 
                    __PACKAGE__ . " has no input with that name, cannot create workflow link!"
            }
        }

        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => $name,
            right_operation => $prediction_operation,
            right_property => $name
        );
    }

    return $prediction_operation;
}

sub _generate_upload_operation {
    my ($self, $workflow, $prediction_operation, %inputs) = @_;

    my $upload_operation;
    if ($self->dump_predictions_to_biosql) {
        $upload_operation = $workflow->add_operation(
            operation_type => Workflow::OperationType::Command->create(
                command_class_name => 'Genome::Model::ProteinAnnotation::Command::UploadResults',
            ),
            name => 'upload results',
        );

        $workflow->add_link(
            left_operation => $prediction_operation,
            left_property => 'bio_seq_features',
            right_operation => $upload_operation,
            right_property => 'bio_seq_features',
        );

        $workflow->add_link(
            left_operation => $upload_operation,
            left_property => 'bio_seq_features',
            right_operation => $workflow->get_output_connector,
            right_property => 'bio_seq_features',
        );
    }
    else {
        $workflow->add_link(
            left_operation => $prediction_operation,
            left_property => 'bio_seq_features',
            right_operation => $workflow->get_output_connector,
            right_property => 'bio_seq_features',
        );
    }

    return $upload_operation;
}

1;

