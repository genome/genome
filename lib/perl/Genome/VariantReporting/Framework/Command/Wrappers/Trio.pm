package Genome::VariantReporting::Framework::Command::Wrappers::Trio;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Framework::Command::Wrappers::Trio {
    is => 'Command::V2',
    has_input => [
        models => {
            is => 'Genome::Model::SomaticValidation',
            is_many => 1,
        },
        snvs_plan_file => {
            is => 'File',
        },
        indels_plan_file => {
            is => 'File',
        },
        output_directory => {
            is => 'Path',
            is_output => 1,
        },
        d0_sample => {
            is => 'Genome::Sample',
            doc => 'Discovery sample',
        },
        d30_sample => {
            is => 'Genome::Sample',
            doc => 'Validation sample',
        },
    ],
};

sub execute {
    my $self = shift;
    my @model_pairs = $self->get_model_pairs;
    for my $model_pair (@model_pairs) {
        my $output_dir = output_dir_for_roi($self->output_directory, $model_pair->{roi});
        my $resource_file = generate_resource_file($model_pair);
        $self->run_reports($model_pair);
    }
    return 1;
}

sub get_model_pairs {
    my $self = shift;
    my %models_for_roi;
    for my $model ($self->models) {
        unless (defined $model->region_of_interest_set_name) {
            $self->warning_message("Skipping model %s because ROI is not defined", $model->__display_name__);
            next;
        }
        push @{$models_for_roi{$model->region_of_interest_set_name}}, $model;
    }

    my @model_pairs;
    for my $roi (keys %models_for_roi) {
        my @models = @{$models_for_roi{$roi}};
        unless (@models == 2) {
            $self->warning_message("Skipping models for ROI %s because there are not exactly two models: %s",
                $roi, join(", ", map {$_->__display_name__} @models));
        }

        my $model_pair;
        $model_pair->{roi} = $roi;

        for my $model (@models) {
            my $build = $model->last_succeeded_build;
            unless (defined $build) {
                $self->warning_message("Skipping model %s because it doesn't have a last succeeded build", $model->__display_name__);
            }
            my $label = $self->model_label($model);
            $model_pair->{$label} = $build;
        }

        unless (defined $model_pair->{discovery} and defined $model_pair->{validation}) {
            $self->warning_message("Skipping model pair for ROI %s because there is not both a discovery and validation model: %s",
                $roi, Data::Dumper::Dumper(\@models));
                next;
        }
        $model_pair->{output_dir} = output_dir_for_roi($self->output_directory, $model_pair->{roi});
        Genome::Sys->create_directory($model_pair->{output_dir});
        push @model_pairs, $model_pair;
    }
    return @model_pairs;
}

sub model_label {
    my $self = shift;
    my $model = shift;
    if ($model->tumor_sample eq $self->d0_sample) {
        return "discovery";
    }
    elsif ($model->tumor_sample eq $self->d30_sample) {
        return "validation";
    }
    else {
        return "unknown";
    }
}

sub output_dir_for_roi {
    my ($output_dir, $roi) = @_;

    return File::Spec->join($output_dir, $roi);
}

sub generate_resource_file {
    my $model_pair = shift;

    my $resource;

    my @aligned_bams;
    push @aligned_bams, $model_pair->{discovery}->merged_alignment_result->id;
    push @aligned_bams, $model_pair->{discovery}->control_merged_alignment_result->id;
    push @aligned_bams, $model_pair->{validation}->control_merged_alignment_result->id;
    $resource->{aligned_bams} = \@aligned_bams;

    my %feature_list_ids;
    my $on_target_feature_list = Genome::FeatureList->get(name => $model_pair->{roi});
    $feature_list_ids{ON_TARGET} = $on_target_feature_list->id;
    my $segdup_feature_list = $model_pair->{discovery}->get_feature_list("segmental_duplications");
    $feature_list_ids{SEG_DUP} = $segdup_feature_list->id;
    $resource->{feature_list_ids} = \%feature_list_ids;

    $resource->{reference_fasta} = $model_pair->{discovery}->reference_sequence_build->full_consensus_path("fa");

    my %translations;
    $translations{d0_tumor} = $model_pair->{discovery}->tumor_sample->name;
    $translations{d30_normal} = $model_pair->{discovery}->normal_sample->name;
    $translations{d30_tumor} = $model_pair->{validation}->tumor_sample->name;
    $resource->{translations} = \%translations;

    YAML::DumpFile(File::Spec->join($model_pair->{output_dir}, "resource.yaml"), $resource);
}

sub run_reports {
    my ($self, $model_pair) = @_;
    Genome::Sys->create_directory($model_pair->{output_dir});

    for my $variant_type(qw(snvs indels)) {
        my $plan_accessor = join('_', $variant_type, "plan_file");
        my %params = (
            input_vcf => resolve_input_vcf($model_pair, $variant_type),
            variant_type => $variant_type,
            output_directory => File::Spec->join($model_pair->{output_dir}, "output"),
            plan_file => $self->$plan_accessor,
            resource_file => File::Spec->join($model_pair->{output_dir}, "resource.yaml"),
            log_directory => File::Spec->join($model_pair->{output_dir}, "logs"),
        );
        print "Params for model pair:\n";
        print Data::Dumper::Dumper(\%params);
        #Genome::VariantReporting::Framework::Command::CreateReport->execute(%params);
    }
}

sub resolve_input_vcf {
    my ($model_pair, $variant_type) = @_;
    $model_pair->{discovery}->get_detailed_vcf_result($variant_type)->get_vcf($variant_type);
}

1;

