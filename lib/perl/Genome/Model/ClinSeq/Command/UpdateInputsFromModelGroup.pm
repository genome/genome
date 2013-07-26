package Genome::Model::ClinSeq::Command::UpdateInputsFromModelGroup;

use strict;
use warnings;

use Genome;
use Command; 
use Data::Dumper;

class Genome::Model::ClinSeq::Command::UpdateInputsFromModelGroup {
    is => 'Command::V2',
    has => [
        update          => { is => 'Genome::Model::ClinSeq',
                            is_many => 1,
                            require_user_verify => 0,
                            doc => 'the models to update' },
        set             => { is => 'Text',
                            is_many => 1,
                            valid_values => ['exome_model'],
                            doc => 'the list of input properties on those models to update' },
        from            => { is => 'Genome::Model', is_many => 1,
                            is_many => 1,
                            require_user_verify => 0,
                            doc => 'the models that should be used as new inputs values' },
        dry_run         => { is => 'Boolean',
                            default_value => 0,
                            is_optional => 1,
                            doc => 'only report possible changes, do not actually update' },
        force           => { is => 'Boolean',
                            default_value => 0,
                            doc => 'allow discrepancies (run first without this flag and be sure you like the changes)' },
    ],
    doc => 'update clinseq models in bulk from another group of models',
};

sub help_synopsis {
    return <<EOS

    # BRCPP models get validation somatic-variation capture data
    genome model clin-seq update-inputs-from-model-group \\
        --update    model_groups.id=59355 \\
        --set       exome_model \\
        --from      model_groups.id=73267 
    
    # HCC models get validation somatic-variation capture data
    genome model clin-seq update-inputs-from-model-group \\
        --update    model_groups.id=73905 \\
        --set       exome_model \\
        --from      model_groups.id=73266

EOS
}

sub execute {
    my $self = shift;
    my @models_to_update = $self->update;
    my @inputs = $self->from;

    my @model_ids_to_update = map { $_->id } @models_to_update;
    
    my %models_to_update_not_matched = map { $_->id => $_ } @models_to_update;
    my %inputs_with_no_matches_wrong_type;
    my %inputs_with_no_matches_correct_type;
    my %inputs_with_multiple_matches;
    my %inputs_matched_with_value_previously_null;
    my %inputs_matched_with_value_previously_different;
    my %inputs_matched_with_value_previously_the_same;

    my %outcomes;
    for my $input (@inputs) {
        my $input_subject = $input->subject;
        my $input_common_name = $input_subject->common_name;
        $self->status_message($input->__display_name__ . ":");

        my %targets;
        if ($input->isa("Genome::Model::SomaticVariation")) {
            my @inst = ($input->tumor_model->instrument_data, $input->normal_model->instrument_data);
            %targets = map { defined($_) ? ($_=>$_) : () } map { $_->target_region_set_name } @inst;
        }

        if ($input->isa("Genome::Model::SomaticVariation") and %targets) {
            # exome/capture somatic model
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
                $self->status_message("\tmultiple candidates!!!!!!!!!!!!!!!!!!!!!!!!");
                $inputs_with_multiple_matches{$input->id} = $input;
                next;
            }
            elsif (@candidates == 0) {
                $self->status_message("\tno candidates!!!!!!!!!!!!!!!!!!!!!!!!");
                $inputs_with_no_matches_correct_type{$input->id} = $input;
                next;
            }
            
            $model_to_update = $candidates[0];    

            $self->status_message("\tfound clinseq model " . $model_to_update->__display_name__);
            if (my $previous_value = $model_to_update->$input_name()) {
                if ($previous_value == $input) {
                    $self->status_message("\tvalue already set!!!!!!!!!!!!!!!!!!!");
                    $inputs_matched_with_value_previously_the_same{$input->id} = $input;
                    next;
                }
                else {
                    $inputs_matched_with_value_previously_different{$input->id} = $input;
                    if(!$self->force) {
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
            $inputs_with_no_matches_wrong_type{$input->id} = $input;
        }
    }

    my %problems = (
        "Models not updated" => \%models_to_update_not_matched,
        "Inputs of the correct type with no matches" => \%inputs_with_no_matches_correct_type,
        "Inputs with multiple matching models" => \%inputs_with_multiple_matches,
        "Inputs matching a model with a value previously set to a different value" => \%inputs_matched_with_value_previously_different,
    );

    my %not_problems = (
        "Inputs matching a model with a value previously the same (ignore)" => \%inputs_matched_with_value_previously_the_same,
        "Inputs matching a model with a value previously null (good)" => \%inputs_matched_with_value_previously_null,
        "Inputs of the wrong type with no matches (ignore)" => \%inputs_with_no_matches_wrong_type,
    );

    my %messages = (%problems, %not_problems);

    $self->status_message("\nSUMMARY:");
    my $errors++;
    for my $msg (reverse sort keys %messages) {
        my $hash = $messages{$msg};
        my @values = values %$hash;
        if ($problems{$msg} and scalar(@values) > 0) {
            $errors++;
            $self->status_message("$msg: " . scalar(@values) . " ********** ");
        }
        else {
            $self->status_message("$msg: " . scalar(@values));
        }
    }

    if ($errors and not $self->force) {
        die $self->error_message("Aborting due to some inconsistencies.  Run with the 'force' option to override.");
    }

    $DB::single = 1;
    return 1;
}

1;

