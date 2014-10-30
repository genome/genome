package Genome::VariantReporting::Command::Wrappers::ModelPairFactory;

use strict;
use warnings;

use Genome;

class Genome::VariantReporting::Command::Wrappers::ModelPairFactory {
    has => {
        models => {
            is => 'Genome::Model::SomaticValidation',
            is_many => 1,
        },
        discovery_sample => { is => 'Genome::Sample', },
        followup_sample => { is => 'Genome::Sample', },
        normal_sample => { is => 'Genome::Sample',},
        output_dir => { is => 'Text', },
        other_input_vcf_pairs => { is => 'Hashref', default_value => {}},
        discovery_output_dir => {
            calculate_from => [qw/output_dir/],
            calculate => q/return File::Spec->join($output_dir, "discovery");/,
        },
        additional_output_dir => {
            calculate_from => [qw/output_dir/],
            calculate => q/return File::Spec->join($output_dir, "followup");/,
        },
        germline_output_dir => {
            calculate_from => [qw/output_dir/],
            calculate => q/return File::Spec->join($output_dir, "germline");/,
        },
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

sub get_model_pairs {
    my $self = shift;

    return if not $self->is_valid;

    my %models_for_roi;
    my @model_pairs;
    for my $model ($self->models) {
        unless (defined $model->region_of_interest_set) {
            $self->warning_message("Skipping model %s because ROI is not defined", $model->__display_name__);
            next;
        }
        if ($self->is_single_bam($model)) {
            push @model_pairs, Genome::VariantReporting::Command::Wrappers::SingleModel->create(
                discovery => $model->last_succeeded_build,
                base_output_dir => $self->germline_output_dir,
            );
        }
        else {
            push @{$models_for_roi{$model->region_of_interest_set->name}}, $model;
        }
    }

    for my $roi (keys %models_for_roi) {

        my @models = @{$models_for_roi{$roi}};
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

        my $discovery_build = $discovery_models[0]->last_succeeded_build;
        if ( not $discovery_build ) {
            $self->warning_message('No last succeeded build for discovery model (%s). Skipping ROI %s.', $discovery_models[0]->__display_name__, $roi);
            next;
        }

        my $validation_build = $validation_models[0]->last_succeeded_build;
        if ( not $validation_build ) {
            $self->warning_message('No last succeeded build for followup model (%s). Skipping ROI %s.', $validation_models[0]->__display_name__, $roi);
            next;
        }

        push @model_pairs, Genome::VariantReporting::Command::Wrappers::ModelPair->create(
            discovery => $discovery_build,
            followup => $validation_build,
            base_output_dir => $self->discovery_output_dir,
        );

        push @model_pairs, Genome::VariantReporting::Command::Wrappers::ModelPair->create(
            discovery => $validation_build,
            followup => $discovery_build,
            base_output_dir => $self->additional_output_dir,
        );

        for my $other_input_vcf_pair (keys %{$self->other_input_vcf_pairs}) {
            push @model_pairs, Genome::VariantReporting::Command::Wrappers::ModelPairWithInput->create(
                discovery => $discovery_build,
                followup => $validation_build,
                plan_file_basename => "cle_docm_report_TYPE.yaml",
                base_output_dir => $self->other_output_dir($other_input_vcf_pair),
                other_snvs_vcf_input => $self->other_input_vcf_pairs->{$other_input_vcf_pair}->[0],
                other_indels_vcf_input => $self->other_input_vcf_pairs->{$other_input_vcf_pair}->[1],
            );
        }
    }

    return @model_pairs;
}

sub other_output_dir {
    my $self = shift;
    my $name = shift;
    return File::Spec->join($self->output_dir, $name);
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

