package Genome::Model::DifferentialExpression;

use strict;
use warnings;

use Genome;

class Genome::Model::DifferentialExpression {
    is => 'Genome::Model',
    has_input => [
        condition_labels_string => {
            is => 'Text',
            doc => 'A comma-delimmited list of labels used to represent each group or condition. ex: "A,B,C,D"',
        },
        condition_model_ids_string => {
            is => 'Text',
            doc => 'Each group is a comma-delimited list of model ids corresponding to a defined label(order matters).  Each group (A, B, C & D) is separated by white space. ex: "1,2 3,4 5,6 7,8"',
        },
    ],
    has_param => [
        transcript_convergence_name => {
            doc => 'The algorithm used to converge transcripts from multiple transcript assemblies(cuffmerge) or the annotation build (cuffcompare).',
            valid_values => ['cuffcompare','cuffmerge'],
            is_optional => 1,
        },
        transcript_convergence_version => {
            doc => 'The version of the transcript convergence algorithm to use.',
            is_optional => 1,
        },
        transcript_convergence_params => {
            doc => 'The parameters used to converge transcripts.',
            is_optional => 1,
        },
        differential_expression_name => {
            doc => 'algorithm used to detect expression levels',
            valid_values => ['cuffdiff'],
            is_optional => 1,
        },
        differential_expression_version => {
            doc => 'the expression detection version used for this model',
            is_optional => 1,
        },
        differential_expression_params => {
            doc => 'the differential expression detection params used for this model',
            is_optional => 1,
        },
        summarize_differential_expression_name => {
            doc => 'algorithm used to summarize differential expression results',
            valid_values => ['cummerbund'],
            is_optional => 1,
        },
        summarize_differential_expression_version => {
            doc => 'The version of the algorithm used to summarize differential expression results',
            is_optional => 1,
        },
    ],
    doc => 'A genome model produced by combining the results of  multiple RNAseq models and determining differential expression of genes between variables.',
};

#TODO: Check the processing_profile, reference_sequence, and eventually annotation to see if they match across models.(Possibly in create?).

#TODO: We need to create a model_group with all the individual models after the model is created but before generating any builds???

#TODO: Resolve the BAM files for each model and map to a string mimicking the condition_model_ids_string

sub _resolve_subject {
    my $self = shift;
    my @subjects = $self->_infer_candidate_subjects_from_input_models();
    if (@subjects > 1) {
        $self->error_message(
            "Conflicting subjects on input models!:\n\t"
            . join("\n\t", map { $_->__display_name__ } @subjects)
        );
        return;
    }
    elsif (@subjects == 0) {
        $self->error_message("No subjects on input models?  Contact Informatics.");
        return;
    }
    return $subjects[0];
    return;
}

sub _infer_candidate_subjects_from_input_models {
    my $self = shift;

    my @condition_model_id_groups = split(/\s+/,$self->condition_model_ids_string);
    my @all_model_ids;
    for my $condition_model_id_group (@condition_model_id_groups) {
        my @condition_model_ids = split(/,/,$condition_model_id_group);
        push @all_model_ids, @condition_model_ids;
    }
    my @input_models;
    for my $model_id (@all_model_ids) {
        my $model = Genome::Model->get($model_id);
        unless ($model) {
            $self->error_message('Failed to get model '. $model_id .' to resolve subject from!');
            return;
        }
        push @input_models, $model;
    }
    my %patients;
    for my $input_model (@input_models) {
        next unless $input_model;
        my $patient;
        if ($input_model->subject->isa("Genome::Individual")) {
            $patient = $input_model->subject;
        }
        else {
            $patient = $input_model->subject->patient;
        }
        $patients{ $patient->id } = $patient;
    }
    my @patients = sort { $a->id cmp $b->id } values %patients;
    if (@patients > 1) {
        $self->status_message('More than one patient found as subject.  Trying the Taxon level.');
        my %taxons;
        for my $patient (@patients) {
            $taxons{ $patient->taxon_id } = $patient->taxon;
        }
        my @taxons = sort { $a->id cmp $b->id } values %taxons;
        return @taxons;
    } else {
        return @patients;
    }
} 
#sub _map_workflow_inputs {
#    my $self = shift;
#    my $build = shift;
#
#    my @inputs = ();
#
#    push @inputs, build_id => $build->id;
#
#    return @inputs;
#}

sub _resolve_workflow_for_build {
    # This is called by Genome::Model::Build::start()
    # Returns a Workflow::Operation
    # By default, builds this from stages(), but can be overridden for custom workflow.
    my $self = shift;
    my $build = shift;
    my $lsf_queue = shift; # TODO: the workflow shouldn't need this yet
    my $lsf_project = shift;

    if (!defined $lsf_queue || $lsf_queue eq '' || $lsf_queue eq 'inline') {
        $lsf_queue = 'apipe';
    }
    if (!defined $lsf_project || $lsf_project eq '') {
        $lsf_project = 'build' . $build->id;
    }

     my @output_properties = qw/
                                  transcript_convergence_result
                              /;
    my $workflow = Workflow::Model->create(
        name => $build->workflow_name,
        input_properties => ['build_id',],
        output_properties => \@output_properties,
    );

    my $log_directory = $build->log_directory;
    $workflow->log_dir($log_directory);

    my $input_connector = $workflow->get_input_connector;
    my $output_connector = $workflow->get_output_connector;

    # Convergence
    my $transcript_convergence_name = $self->transcript_convergence_name;
    my $transcript_convergence_operation;
    if ($self->transcript_convergence_name eq 'cuffcompare') {
        $transcript_convergence_operation = $workflow->add_operation(
            name => 'Differential Expression',
            operation_type => Workflow::OperationType::Command->create(
                command_class_name => 'Genome::Model::DifferentialExpression::Command::Cuffcompare',
            )
        );
    } elsif ($self->transcript_convergence_name eq 'cuffmerge') {
        $transcript_convergence_operation = $workflow->add_operation(
            name => 'Differential Expression',
            operation_type => Workflow::OperationType::Command->create(
                command_class_name => 'Genome::Model::DifferentialExpression::Command::Cuffmerge',
            )
        );
    } else {
        die('Unsupported transcript_convergence_name: '. $self->transcript_convergence_name);
    }
    
    $transcript_convergence_operation->operation_type->lsf_queue($lsf_queue);
    $transcript_convergence_operation->operation_type->lsf_project($lsf_project);
    
    # Differential Expression

    # Summary Results

    return $workflow;
}


1;

