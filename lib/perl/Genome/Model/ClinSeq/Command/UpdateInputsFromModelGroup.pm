package Genome::Model::ClinSeq::Command::UpdateInputsFromModelGroup;

use strict;
use warnings;

use Genome;
use Command; 
use Data::Dumper;

class Genome::Model::ClinSeq::Command::UpdateInputsFromModelGroup {
    is => 'Command::V2',
    has => [
        clinseq_models  => { is => 'Genome::Model::ClinSeq',
                            is_many => 1,
                            shell_args_position => 1,
                            require_user_verify => 0,
                            doc => 'the clinseq models to update' },
        input_models    => { is => 'Genome::Model', is_many => 1,
                            is_many => 1,
                            shell_args_position => 2,
                            require_user_verify => 0,
                            doc => 'the models that should be used as inputs to the clinseq models' },
        set             => { is => 'Text',
                            is_many => 1,
                            is_optional => 1,
                            valid_values => ['exome_model'],
                            doc => 'the list of properties to update on the clinseq model, by default any for which the data type fits' },
        dry_run         => { is => 'Boolean',
                            default_value => 0,
                            is_optional => 1,
                            doc => 'only report possible changes, do not actually update' },
        allow_previous  => { is => 'Boolean',
                            default_value => 0,
                            doc => 'set the input even if it is currently set (by default only null inputs are set)' },
        allow_unmatched  => { is => 'Boolean',
                            default_value => 0,
                            doc => 'allow some input models that do not match the clinseq models, or clinseq models that do not match input models' },
    ],
    doc => 'update clinseq models in bulk from another group of models',
};

sub help_synopsis {
    return <<EOS
    genome model clin-seq update-inputs-from-model-group \
        --clinseq-models    model_groups.name=myclinseqgroup \
        --input-models      model_groups.name=myexomegroup \
    
    genome model clin-seq update-inputs-from-model-group \
        --clinseq-models    model_groups.name=myclinseqgroup \
        --input-models      model_groups.name=myrnaseqgroup \
EOS
}

sub execute {
    my $self = shift;
    my @models_to_update = $self->clinseq_models;
    my @inputs = $self->input_models;

    my @model_ids_to_update = map { $_->id } @models_to_update;
    
    my %models_to_update_not_matched = map { $_->id => $_ } @models_to_update;
    my %inputs_with_no_matches;
    my %inputs_with_multiple_matches;
    my %inputs_matched_with_value_previously_null;
    my %inputs_matched_with_value_previously_different;
    my %inputs_matched_with_value_previously_the_same;

    my %outcomes;
    for my $input (@inputs) {
        my $input_subject = $input->subject;
        my $input_common_name = $input_subject->common_name;
        $self->status_message($input->__display_name__ . ":");

        if ($input->isa("Genome::Model::SomaticVariation") and $input->tumor_model->target_region_set_name) {
            # exome/capture
            my $model_to_update;
            my $input_name = 'exome_model';

            # approach 1: find the clinseq model to update by finding WGS models with the same subject
            my @candidates = Genome::Model::ClinSeq->get(
                id => \@model_ids_to_update,
                "wgs_model.subject.id" => $input_subject->id 
            );
            for my $candidate (@candidates) {
                delete $models_to_update_not_matched{$candidate->id};
            }
            if (@candidates > 1) {
                $self->warning_message("\tmultiple candidates!");
                $inputs_with_multiple_matches{$input->id} = $input;
                next;
            }
            elsif (@candidates == 0) {
                $self->warning_message("\tno candidates!");
                $inputs_with_no_matches{$input->id} = $input;
                next;
            }
            
            $model_to_update = $candidates[0];    

            $self->status_message("\tfound clinseq model " . $model_to_update->__display_name__);
            if (my $previous_value = $model_to_update->$input_name()) {
                if ($previous_value == $input) {
                    $self->status_message("\tvalue already set!");
                    $inputs_matched_with_value_previously_the_same{$input->id} = $input;
                    next;
                }
                else {
                    $inputs_matched_with_value_previously_different{$input->id} = $input;
                    if(!$self->allow_previous) {
                        $self->status_message("\tSKIP UPDATE OLD $input_name on " . $model_to_update->__display_name__ . " from " . $previous_value->__display_name__ . " to " . $input->__display_name__ . ".");
                        next;
                    }
                    else {
                        $self->status_message("\tUPDATE OLD $input_name on " . $model_to_update->__display_name__ . " from " . $previous_value->__display_name__ . " to " . $input->__display_name__ . ".");
                    }
                }
            }
            else {
                $inputs_matched_with_value_previously_null{$input->id} = $input;
                $self->status_message("\tset $input_name on " . $model_to_update->__display_name__ . " to " . $input->__display_name__ . ".");
            }

            unless ($self->dry_run) {
                $model_to_update->$input_name($input);
            }
        }
        else {
            die "support only exists for somatic variation models.  Add support for more types here: " . __FILE__ . ", line " . __LINE__;
        }
    }

    my %problems = (
        "Models not updated" => \%models_to_update_not_matched,
        "Inputs with no matching models" => \%inputs_with_no_matches,
        "Inputs with multiple matching models" => \%inputs_with_multiple_matches,
        "Inputs matching a model with a value previously set to a different value" => \%inputs_matched_with_value_previously_different,
    );

    my %not_problems = (
        "Inputs matching a model with a value previously the same" => \%inputs_matched_with_value_previously_the_same,
        "Inputs matching a model with a value previously null" => \%inputs_matched_with_value_previously_null,
    );

    my %messages = (%problems, %not_problems);

    my $errors++;
    for my $msg (reverse sort keys %messages) {
        my $hash = $messages{$msg};
        my @values = values %$hash;
        $self->status_message("$msg: " . scalar(@values));
        if ($problems{$msg}) {
            $errors++;
        }
    }

    if ($errors) {
        die $self->error_message("Aborting due to some inconsistencies.  Run with the 'force' option to override.");
    }

    return 1;
}

1;

