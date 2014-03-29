package Genome::Model::Command::TumorNormalDefine;

use strict;
use warnings;

use Genome;
use Memoize;

class Genome::Model::Command::TumorNormalDefine {
    is => 'Genome::Command::Base',
    has_input => [
        reference_alignment_models => {
            is => 'Genome::Model',
            is_many => 1,
            shell_args_position => 1,
            require_user_verify => 0,
            doc => 'A list of reference-alignment models or a model-group',
        },
        processing_profile => {
            is => 'Genome::ProcessingProfile',
            doc => 'The processing-profile that will be used when defining '.
                'the new model(s)',
        },
        annotation_build => {
            is => 'Genome::Model::Build::ImportedAnnotation',
            doc => 'The annotation-build that will be used when defining the '.
                'new model(s)',
        },
        model_group_name => {
            is => 'Text',
            doc => 'The name of the model group that will be created',
        },
    ],
    has_optional_input => [
        previously_discovered_variations_build => {
            is => 'Genome::Model::Build::ImportedVariationList',
            doc => 'A build with previously discovered variants that will '.
                'be used when defining the new model(s)',
        },
        identify => {
            is => 'Text',
            valid_values => ['tumor', 'normal'],
            default => 'normal',
            doc => 'Determines which (tumor or normal) group is selected by '.
                    'matching the regex',
        },
        regex => {
            is => 'Text',
            doc => "A regular expression that will match against the ".
                    "extraction_label of the models' subject.  If not ".
                    "specified, the user will be interactively prompted."
        },
        interactive => {
            is => 'Boolean',
            default => 1,
            doc => 'Interactively prompt the user.',
        },
    ],
    doc => 'Define new somatic-variation model(s) from a set of '.
        'reference-alignment models.'
};

sub help_synopsis {
    return <<EOS

Define new somatic-variation model(s) from a set of reference-alignment models.
EOS
}

sub help_detail {
    return <<EOS

The <regex> input is used to determine the models which are either the tumor
models or the normal models.  For example with this set of models:

H_NI-LBHDec-014-LBHDec-014_skin.prod-refalign.capture.AML
H_NI-LBHDec-014-LBHDec-014_C111.prod-refalign.capture.AML
H_NI-LBHDec-014-LBHDec-014_C222.prod-refalign.capture.AML
H_NI-LBHDec-014-LBHDec-014_C333.prod-refalign.capture.AML

One might define the normal model with the regex 'skin'.  Then the tumor models
would be the last three.  This command would then create 3 somatic-variation
models and put them into a model group.

Alternatively, with the same set of models as inputs one might define the tumor
models with the regex 'C[23]*'.  Then the last two models would be identified
as tumor and the first two identified as normal models.  This command would
then create 4 somatic-variation models (one for every tumor-normal combination)
and put them into a model group.
EOS
}

sub execute {
    my $self = shift;

    my $model_group = Genome::ModelGroup->create(
        name => $self->model_group_name,
    );

    unless ($model_group) {
        $self->error_message('Failed to create model group for name: '.
            $self->model_group_name);
        return;
    }
    $self->status_message('Created model group:');
    $self->status_message('ID: ' . $model_group->id . ', NAME: ' . $model_group->name);

    my @reference_alignment_models = $self->reference_alignment_models;
    $self->status_message("Found %d Reference Alignment Models.",
        scalar(@reference_alignment_models));

    $self->resolve_regular_expression();

    $self->status_message("This will define %d new SomaticVariation models",
        scalar($self->define_params));
    if ($self->interactive) {
        return unless prompt_yn("Are you sure you want to continue?");
    }

    my @new_models;
    for my $params ($self->define_params) {
        my $cmd = Genome::Model::Command::Define::SomaticVariation->create(
            %{$params});
        unless ($cmd->execute()) {
            die "Couldn't define SomaticVariation model with params " .
                Data::Dumper::Dumper($params);
        }
        my $new_model = Genome::Model->get(name => $params->{model_name});
        push @new_models, $new_model;
    }

    $self->status_message('Assigning %d to group (%s).',
        scalar(@new_models),
        $model_group->name,
    );
    $model_group->assign_models(@new_models);

    return 1;
}

sub linked_models {
    my $self = shift;

    my ($match_label, $non_match_label);
    if ($self->identify eq 'tumor') {
        $match_label = 'tumor';
        $non_match_label = 'normal';
    } else {
        $match_label = 'normal';
        $non_match_label = 'tumor';
    }

    my %result;
    for my $model ($self->matching_models) {
        my $patient_name = $model->subject->patient_name;

        $result{$patient_name}{$match_label} =
            [grep {$_->subject->patient_name eq $patient_name} $self->matching_models];
        $result{$patient_name}{$non_match_label} =
            [grep {$_->subject->patient_name eq $patient_name} $self->non_matching_models];
    }
    return %result;
}

sub matching_models {
    my $self = shift;

    my $regex = $self->regex;
    return grep {$_->subject->extraction_label =~ /$regex/} $self->reference_alignment_models;
}
Memoize::memoize('matching_models');

sub non_matching_models {
    my $self = shift;

    my $regex = $self->regex;
    return grep {not $_->subject->extraction_label =~ /$regex/} $self->reference_alignment_models;
}
Memoize::memoize('non_matching_models');

sub base_define_params {
    my $self = shift;

    my %params = (
        processing_profile => $self->processing_profile,
        annotation_build => $self->annotation_build,
    );
    if ($self->previously_discovered_variations_build) {
        $params{previously_discovered_variations_build} =
            $self->previously_discovered_variations_build;
    }
    return %params;
}

sub define_params {
    my $self = shift;
    my %linked_models = $self->linked_models;

    my @result;
    for my $models (values %linked_models) {
        for my $tumor_model (@{$models->{'tumor'}}) {
            for my $normal_model (@{$models->{'normal'}}) {
                my $model_name = sprintf("%s (%s vs %s)",
                    $tumor_model->subject->patient_name,
                    $tumor_model->subject->extraction_label,
                    $normal_model->subject->extraction_label,
                );
                push @result, {
                    'tumor_model' => $tumor_model,
                    'normal_model' => $normal_model,
                    'model_name' => $model_name,
                    $self->base_define_params,
                };
            }
        }
    }
    return @result;
}
Memoize::memoize('define_params');

sub resolve_regular_expression {
    my $self = shift;

    if ($self->regex and !$self->interactive) {
        $self->display_matches();
    } else {
        my $done;
        if ($self->regex) {
            $done = $self->confirm_regex();
        }

        while (!$done) {
            $self->regex(
                prompt(sprintf("Enter regex that identifies a %s: ", $self->identify))
            );
            $done = $self->confirm_regex();
        }
    }
}

sub confirm_regex {
    my $self = shift;

    $self->display_matches();
    return prompt_yn(sprintf("Use regex /%s/ to identify %s?",
        $self->regex, $self->identify));
}

sub display_matches {
    my $self = shift;

    my $regex = $self->regex;

    my @matches;
    my @misses;
    for my $label ($self->extraction_labels) {
        if ($label =~ m/$regex/) {
            push @matches, $label;
        } else {
            push @misses, $label;
        }
    }

    $self->status_message('Regex /%s/ matches %d of %d extraction labels.',
        $regex, scalar(@matches), scalar($self->extraction_labels));

    if (scalar(@matches)) {
        $self->status_message('  Example Match: %s', $matches[0]);
        if (scalar(@misses)) {
            $self->status_message('  Example Miss: %s', $misses[0]);
        }
    }
}

sub prompt_yn {
    my $query = shift;
    my $answer = prompt("$query (y/N): ");
    return lc($answer) eq 'y';
}

sub prompt {
    my $query = shift;
    local $| = 1;
    print $query;
    chomp(my $answer = <STDIN>);
    return $answer;
}

sub extraction_labels {
    my $self = shift;
    return map {$_->subject->extraction_label} $self->reference_alignment_models;
}
Memoize::memoize('extraction_labels');

1;
