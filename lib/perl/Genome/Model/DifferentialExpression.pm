package Genome::Model::DifferentialExpression;

use strict;
use warnings;

use Genome;

class Genome::Model::DifferentialExpression {
    is => 'Genome::Model',
    has_input => [
        condition_pairs => {
            is => 'Text',
            is_many => 1,
            doc => 'Each condition-pair is a space-delimited pair consisting of the label and a model-id. Labels should not contain spaces. Example: "A 1","A 2","B 3" will associate label A with models 1 and 2, and label B with model 3',
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

# Creates a sorted mapping of labels to model ids
sub condition_pairs_sorted_mapping {
    my ($self) = @_;

    my %hash = $self->condition_pairs_unsorted_hash();

    my @mapping; # sort order should be consistent
    my @sorted_labels = sort keys %hash;
    for my $label (@sorted_labels) {
        my @sorted_ids = sort @{$hash{$label}};
        push @mapping, {
            label     => $label,
            model_ids => [@sorted_ids],
        };
    }

    return @mapping;
}

sub condition_pairs_unsorted_hash {
    my ($self) = @_;

    my %hash;
    my @pairs = $self->condition_pairs();
    for my $pair (@pairs) {
        my $regex = qr/^
          \s*             # strip leading whitespace
          (\S+            # match at least one word as the label
            (?:\s+\S+)*   # (optionally several more words)
          )
          \s+             # match whitespace
          ([0-9a-fA-F]+)  # match an ID
          \s*             # strip trailing whitespace
        $/x;
        unless ($pair =~ $regex) {
            die $self->error_message(
                "Could not extract a model-id and label from " .
                "'$pair' (expected format: 'label model-id')"
            );
        }
        my ($label, $model_id) = ($1, $2);

        push @{$hash{$label}}, $model_id;
    }

    return %hash;
}

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
            $patient = $input_model->subject->individual;
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

    my @pairs_mapping = $self->condition_pairs_sorted_mapping;
    my @all_model_ids = map { @{$_->{model_ids}} } @pairs_mapping;
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
    my $self = shift;
    my $build = shift;
    my $lsf_queue = shift; # TODO: the workflow shouldn't need this yet
    my $lsf_project = shift;

    if (!defined $lsf_queue || $lsf_queue eq '' || $lsf_queue eq 'inline') {
        $lsf_queue = Genome::Config::get('lsf_queue_build_worker_alt');
    }
    if (!defined $lsf_project || $lsf_project eq '') {
        $lsf_project = 'build' . $build->id;
    }

    my $workflow = Genome::WorkflowBuilder::DAG->create(
        name => $build->workflow_name,
    );

    #TODO: Validate that all models have a Succeeded build or at minimum an Alignment Result and Cufflinks output
    unless ($self->_validate_rna_seq_succeeded_build) {
        die "One or more input rnaseq models do not have a succeeded build. "
            ."Please check your input models to make sure they have all succeeded.";
    }

    # Convergence
    my $transcript_convergence_name = $self->transcript_convergence_name;
    my $transcript_convergence_operation;
    if ($self->transcript_convergence_name eq 'cuffcompare') {
        $transcript_convergence_operation = Genome::WorkflowBuilder::Command->create(
            name => 'Transcript Convergence',
            command => 'Genome::Model::DifferentialExpression::Command::Cuffcompare',
        );
        $workflow->add_operation($transcript_convergence_operation);
    } elsif ($self->transcript_convergence_name eq 'cuffmerge') {
        $transcript_convergence_operation = Genome::WorkflowBuilder::Command->create(
            name => 'Transcript Convergence',
            command => 'Genome::Model::DifferentialExpression::Command::Cuffmerge',
        );
        $workflow->add_operation($transcript_convergence_operation);
    } else {
        die('Unsupported transcript_convergence_name: '. $self->transcript_convergence_name);
    }

    $transcript_convergence_operation->lsf_queue($lsf_queue);
    $transcript_convergence_operation->lsf_project($lsf_project);
    $workflow->connect_input(
        input_property => 'build_id',
        destination => $transcript_convergence_operation,
        destination_property => 'build_id'
    );

    # Differential Expression
    my $differential_expression_name = $self->differential_expression_name;
    my $differential_expression_operation;
    if ($self->differential_expression_name eq 'cuffdiff') {
        $differential_expression_operation = Genome::WorkflowBuilder::Command->create(
            name => 'Differential Expression',
            command => 'Genome::Model::DifferentialExpression::Command::Cuffdiff',
        );
        $workflow->add_operation($differential_expression_operation);
     } else {
         die('Unsupported differential_expression_name: '. $self->differential_expression_name);
    }

    $differential_expression_operation->lsf_queue($lsf_queue);
    $differential_expression_operation->lsf_project($lsf_project);
    $workflow->create_link(
        source => $transcript_convergence_operation,
        source_property => 'build_id',
        destination => $differential_expression_operation,
        destination_property => 'build_id'
    );

    # Summarize Results
    my $summarize_name = $self->summarize_differential_expression_name;
    if ($summarize_name) {
        my $summarize_operation;
        if ($summarize_name eq 'cummerbund') {
            $summarize_operation = Genome::WorkflowBuilder::DAG->create(
                name => 'Summarize Differential Expression',
                command => 'Genome::Model::DifferentialExpression::Command::Cummerbund',
            );
            $workflow->add_operation($summarize_operation);
        } else {
            die('Unsupported summarize differential_expression_name: '. $summarize_name);
        }
        $summarize_operation->lsf_queue($lsf_queue);
        $summarize_operation->lsf_project($lsf_project);
        $workflow->create_link(
            source => $differential_expression_operation,
            source_property => 'build_id',
            destination => $summarize_operation,
            destination_property => 'build_id'
        );
        $workflow->connect_output(
            source => $summarize_operation,
            source_property => 'result',
            output_property => 'summarize_result'
        );
    }

    $workflow->connect_output(
        source => $transcript_convergence_operation,
        source_property => 'result',
        output_property => 'transcript_convergence_result'
    );
    $workflow->connect_output(
        source => $differential_expression_operation,
        source_property => 'result',
        output_property => 'differential_expression_result'
    );

    my $log_directory = $build->log_directory;
    $workflow->recursively_set_log_dir($log_directory);

    return $workflow;
}


1;

