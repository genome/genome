package Genome::VariantReporting::Command::Wrappers::ModelPairFactory;

use strict;
use warnings;

use Genome;

class Genome::VariantReporting::Command::Wrappers::ModelPairFactory {
    is => ['Genome::VariantReporting::Command::Wrappers::Utils'],
    has => {
        models => {
            is => 'Genome::Model::SomaticValidation',
            is_many => 1,
        },
        discovery_sample => { is => 'Genome::Sample', },
        followup_sample => { is => 'Genome::Sample', is_optional => 1},
        normal_sample => { is => 'Genome::Sample',},
        other_input_vcf_pairs => { is => 'Hashref', default_value => {}},
    },
};

sub is_valid {
    my $self = shift;

    if (my @problems = $self->__errors__) {
        $self->error_message('Model pair factory is invalid!');
        for my $problem (@problems) {
            my @properties = $problem->properties;
            $self->error_message("Property " .
                join(',', map { "'$_'" } @properties) .
                ': ' . $problem->desc);
        }
        return;
    }

    return 1;
}

sub get_models_for_roi {
    my $self = shift;

    my %models_for_roi;
    for my $model ($self->models) {
        next unless defined $model->region_of_interest_set;
        push @{$models_for_roi{$model->region_of_interest_set->name}}, $model;
    }
    return \%models_for_roi;
}

sub get_model_pairs {
    my $self = shift;

    return if not $self->is_valid;

    my @model_pairs;
    for my $model_list (values %{$self->get_models_for_roi}) {
        for my $model (@{$model_list}) {
            if ($self->is_single_bam($model)) {
                push @model_pairs, Genome::VariantReporting::Command::Wrappers::SingleModel->create(
                    common_translations => $self->get_common_translations(),
                    discovery => $model->last_succeeded_build,
                    label => 'germline',
                );
            }
        }
    }

    my %models_for_roi = %{$self->get_models_for_roi};
    while (my ($roi, $model_list) = each %models_for_roi) {
        my @models = grep {!$self->is_single_bam($_)} @{$model_list};
        my ($discovery_build, $validation_build);

        if (@models == 1) {#case for no followup
            $discovery_build = $models[0]->last_succeeded_build;
            push @model_pairs, Genome::VariantReporting::Command::Wrappers::ModelPair->create(
                common_translations => $self->get_common_translations(),
                plan_file_basename => "cle_somatic_TYPE_report.yaml",
                discovery => $discovery_build,
                label => 'discovery',
            );
        }
        else {
            unless (@models == 2) {
                $self->warning_message("Skipping models for ROI %s because there are not exactly two models: %s",
                    $roi, join(", ", map {$_->__display_name__} @models));
                next;
            }

            my @discovery_models = grep { $self->is_model_discovery($_) } @models;
            my @validation_models = grep { $self->is_model_followup($_) } @models;

            if ( @discovery_models != 1 or @validation_models != 1 ) {
                $self->warning_message("Incorrect discovery/followup pairing for models for ROI (%s). One of each is required!\nDiscovery:%s\nFollowup:%s\n", $roi, join(", ", map {$_->__display_name__} @discovery_models),
                    join(", ", map {$_->__display_name__} @validation_models));
                next;
            }

            $discovery_build = $discovery_models[0]->last_succeeded_build;
            if ( not $discovery_build ) {
                $self->warning_message('No last succeeded build for discovery model (%s). Skipping ROI %s.', $discovery_models[0]->__display_name__, $roi);
                next;
            }

            $validation_build = $validation_models[0]->last_succeeded_build;
            if ( not $validation_build ) {
                $self->warning_message('No last succeeded build for followup model (%s). Skipping ROI %s.', $validation_models[0]->__display_name__, $roi);
                next;
            }

            push @model_pairs, Genome::VariantReporting::Command::Wrappers::ModelPair->create(
                common_translations => $self->get_common_translations(),
                discovery => $discovery_build,
                followup => $validation_build,
                label => "discovery",
            );

            push @model_pairs, Genome::VariantReporting::Command::Wrappers::ModelPair->create(
                common_translations => $self->get_common_translations(),
                discovery => $validation_build,
                followup => $discovery_build,
                label => "followup",
            );
        }

        for my $other_input_vcf_pair (keys %{$self->other_input_vcf_pairs}) {
            my %params = (
                common_translations => $self->get_common_translations(),
                discovery => $discovery_build,
                plan_file_basename => "cle_docm_report_TYPE.yaml",
                label => $other_input_vcf_pair,
                other_snvs_vcf_input => $self->other_input_vcf_pairs->{$other_input_vcf_pair}->[0],
                other_indels_vcf_input => $self->other_input_vcf_pairs->{$other_input_vcf_pair}->[1],
            );
            if ($validation_build) {
                $params{followup} = $validation_build
            }
            else { # case for no followup
                $params{plan_file_basename} = "cle_somatic_docm_report_TYPE.yaml";
            }
            push @model_pairs, Genome::VariantReporting::Command::Wrappers::ModelPairWithInput->create(%params);
        }
    }

    return \@model_pairs;
}

sub get_common_translations {
    my $self = shift;

    my $sample_translations  = {
        $self->discovery_sample->name => sprintf('Discovery(%s)', $self->discovery_sample->name),
        $self->normal_sample->name    => sprintf('Normal(%s)', $self->normal_sample->name),
    };

    my $library_translations = {
        $self->get_library_name_labels('discovery', $self->discovery_sample, [$self->models]),
        $self->get_library_name_labels('normal', $self->normal_sample, [$self->models]),
    };

    if ($self->followup_sample) {
        $sample_translations->{$self->followup_sample->name} = sprintf('Followup(%s)', $self->followup_sample->name);
        %$library_translations = (%$library_translations, $self->get_library_name_labels('followup', $self->followup_sample, [$self->models]));
    }
    return {
        sample_name_labels  => $sample_translations,
        library_name_labels => $library_translations,
    };
}

sub is_model_discovery {
    my ($self, $model) = @_;
    return $self->discovery_sample->id eq $model->tumor_sample->id;
}

sub is_model_followup {
    my ($self, $model) = @_;
    return $self->followup_sample->id eq $model->tumor_sample->id;
}

sub is_single_bam {
    my ($self, $model) = @_;
    return (!(defined $model->normal_sample) and $self->normal_sample->id eq $model->tumor_sample->id);
}

1;
