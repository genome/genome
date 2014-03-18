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
        reference_sequence_build => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            is_optional => 1,
        },
        annotation_build => {
            is => "Genome::Model::Build::ImportedAnnotation",
            is_optional => 1,
        }
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
        transcript_convergence_biotypes => {
            doc => 'The transcript feature biotypes to include in the convergence when using cuffcompare. example: protein_coding,pseudogene,miRNA,lincRNA,snoRNA,snRNA',
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
        differential_expression_mask_reference_transcripts => {
            doc => 'the annotation file basename used to mask known reference transcripts during differential expression.',
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

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    unless ($self) { return; }

    unless ($self->_validate_rna_seq_processing_profile) {
        return;
    }
    unless ($self->reference_sequence_build) {
        unless ($self->_resolve_reference_sequence_build) {
            return;
        }
    }
    unless ($self->annotation_build) {
        unless ($self->_resolve_annotation_build) {
            return;
        }
    }
    return $self;
}

sub _validate_rna_seq_succeeded_build {
    my $self = shift;

    for my $input_model ($self->input_models) {
        unless ($input_model->last_succeeded_build) {
            $self->error_message('The input model '. $input_model->id. ' does not have a succeeded build!');
            return;
        }
    }
    
    return 1;
}

sub _validate_rna_seq_processing_profile {
    my $self = shift;
    
    my %processing_profiles = map {$_->processing_profile->id => $_->processing_profile } $self->input_models;
    my @pp_ids = keys %processing_profiles;
    if (scalar(@pp_ids) == 1) {
        return 1;
    }
    $self->error_message('Multiple processing profile ids found for input models: '. "\n". join("\n",@pp_ids));
    return;
}

sub _resolve_reference_sequence_build {
    my $self = shift;
    
    my %reference_sequence_builds = map {$_->reference_sequence_build->id => $_->reference_sequence_build} $self->input_models;
    my @build_ids = keys %reference_sequence_builds;
    if (scalar(@build_ids) == 1) {
        $self->reference_sequence_build($reference_sequence_builds{$build_ids[0]});
        return 1;
    }
    $self->error_message('Multiple reference sequence build ids found for input models: '. "\n". join("\n",@build_ids));
    return;
}

sub _resolve_annotation_build {
    my $self = shift;
    
    my %annotation_builds = map {$_->annotation_build->id => $_->annotation_build} $self->input_models;
    my @build_ids = keys %annotation_builds;
    if (scalar(@build_ids) == 1) {
        $self->annotation_build($annotation_builds{$build_ids[0]});
        return 1;
    }
    $self->error_message('Multiple annotation build ids found for input models: '. "\n". join("\n",@build_ids));
    return;

}

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
    my $subject = $subjects[0];
    return $subject;
}

sub _infer_candidate_subjects_from_input_models {
    my $self = shift;

    my %patients;
    for my $input_model ($self->input_models) {
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

sub input_models {
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
            $self->error_message('Failed to get model '. $model_id .' as input model to '. __PACKAGE__);
            return;
        }
        push @input_models, $model;
    }
    return @input_models;
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
        $lsf_queue = $ENV{GENOME_LSF_QUEUE_BUILD_WORKER_ALT};
    }
    if (!defined $lsf_project || $lsf_project eq '') {
        $lsf_project = 'build' . $build->id;
    }

     my @output_properties = qw/
                                  transcript_convergence_result
                                  differential_expression_result
                                  summarize_result
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

    #TODO: Validate that all models have a Succeeded build or at minimum an Alignment Result and Cufflinks output
    unless ($self->_validate_rna_seq_succeeded_build) {
        die "One or more input rnaseq models do not have a succeeded build. "
            ."Please check your input models to make sure they have all succeeded.";
    }

    # Convergence
    my $transcript_convergence_name = $self->transcript_convergence_name;
    my $transcript_convergence_operation;
    if ($self->transcript_convergence_name eq 'cuffcompare') {
        $transcript_convergence_operation = $workflow->add_operation(
            name => 'Transcript Convergence',
            operation_type => Workflow::OperationType::Command->create(
                command_class_name => 'Genome::Model::DifferentialExpression::Command::Cuffcompare',
            )
        );
    } elsif ($self->transcript_convergence_name eq 'cuffmerge') {
        $transcript_convergence_operation = $workflow->add_operation(
            name => 'Transcript Convergence',
            operation_type => Workflow::OperationType::Command->create(
                command_class_name => 'Genome::Model::DifferentialExpression::Command::Cuffmerge',
            )
        );
    } else {
        die('Unsupported transcript_convergence_name: '. $self->transcript_convergence_name);
    }
    
    $transcript_convergence_operation->operation_type->lsf_queue($lsf_queue);
    $transcript_convergence_operation->operation_type->lsf_project($lsf_project);
    $workflow->add_link(
        left_operation => $input_connector,
        left_property => 'build_id',
        right_operation => $transcript_convergence_operation,
        right_property => 'build_id'
    );
    
    # Differential Expression
    my $differential_expression_name = $self->differential_expression_name;
    my $differential_expression_operation;
    if ($self->differential_expression_name eq 'cuffdiff') {
        $differential_expression_operation = $workflow->add_operation(
            name => 'Differential Expression',
            operation_type => Workflow::OperationType::Command->create(
                command_class_name => 'Genome::Model::DifferentialExpression::Command::Cuffdiff',
            )
        );
     } else {
         die('Unsupported differential_expression_name: '. $self->differential_expression_name);
    }
    
    $differential_expression_operation->operation_type->lsf_queue($lsf_queue);
    $differential_expression_operation->operation_type->lsf_project($lsf_project);
    $workflow->add_link(
        left_operation => $transcript_convergence_operation,
        left_property => 'build_id',
        right_operation => $differential_expression_operation,
        right_property => 'build_id'
    );
    
    # Summarize Results
    # TODO: cummerbund
    my $summarize_name = $self->summarize_differential_expression_name;
    my $summarize_operation;
    if ($summarize_name eq 'cummerbund') {
        $summarize_operation = $workflow->add_operation(
            name => 'Summarize Differential Expression',
            operation_type => Workflow::OperationType::Command->create(
                command_class_name => 'Genome::Model::DifferentialExpression::Command::Cummerbund',
            )
        );
    } else {
         die('Unsupported summarize differential_expression_name: '. $summarize_name);
     }
    
    $summarize_operation->operation_type->lsf_queue($lsf_queue);
    $summarize_operation->operation_type->lsf_project($lsf_project);
    $workflow->add_link(
        left_operation => $differential_expression_operation,
        left_property => 'build_id',
        right_operation => $summarize_operation,
        right_property => 'build_id'
    );

    $workflow->add_link(
        left_operation => $summarize_operation,
        left_property => 'result',
        right_operation => $output_connector,
        right_property => 'summarize_result'
    );
    $workflow->add_link(
        left_operation => $transcript_convergence_operation,
        left_property => 'result',
        right_operation => $output_connector,
        right_property => 'transcript_convergence_result'
    );
    $workflow->add_link(
        left_operation => $differential_expression_operation,
        left_property => 'result',
        right_operation => $output_connector,
        right_property => 'differential_expression_result'
    );
    return $workflow;
}


1;

