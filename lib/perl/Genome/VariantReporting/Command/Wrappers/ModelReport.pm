package Genome::VariantReporting::Command::Wrappers::ModelReport;

use strict;
use warnings;

use Genome;

use File::Basename qw(dirname);
use File::Spec;

class Genome::VariantReporting::Command::Wrappers::ModelReport {
    is => ['Command::V2', 'Genome::VariantReporting::Command::Wrappers::Utils'],
    has_input => {
        model => {
            is => 'Genome::Model::SomaticValidation',
        },
    },

    has => [
        normal_sample => {
          is => 'Genome::Sample',
          via => 'model',
        },

        tumor_sample => {
          is => 'Genome::Sample',
          via => 'model',
        },
    ],
};

sub execute {
    my $self = shift;

    my $model = $self->model;
     my $model_pair;
    if ($self->is_single_bam($model)) {
        # Germline
        $model_pair = Genome::VariantReporting::Command::Wrappers::SingleModel->create(
            common_translations => $self->get_germline_translations,
            discovery => $self->build,
            label => 'germline',
        );
    } else {
        #Somatic
        $model_pair = Genome::VariantReporting::Command::Wrappers::ModelPair->create(
            common_translations => $self->get_somatic_translations,
            discovery => $self->build,
            plan_file_basename => 'somatic_TYPE_report.yaml',
            label => 'somatic',
        );
    }
    for my $variant_type(qw(snvs indels)) {
        my %params = (
            input_vcf => $model_pair->input_vcf($variant_type),
            variant_type => $variant_type,
            plan_file => $model_pair->plan_file($variant_type),
            translations_file => $model_pair->translations_file,
        );
        my $cmd = Genome::VariantReporting::Command::CreateReport->create(%params);

        my $process = $cmd->execute();
        $process->add_note(
            header_text => 'creation metadata',
            body_text => sprintf('%s report created by "genome variant-reporting wrappers model-report --model=%s"',
                $variant_type, $self->model->id),
        );
    }
    return 1;
};

sub get_germline_translations {
    my $self = shift;

    return {
        sample_name_labels => {
            $self->get_sample_name_labels('tumor'),
        },
        library_name_labels => {
            $self->get_library_name_labels('tumor', $self->tumor_sample, [$self->build]),
        },
    };
}

sub get_somatic_translations {
    my $self = shift;

    return {
        sample_name_labels => {
            $self->get_sample_name_labels('tumor'),
            $self->get_sample_name_labels('normal'),
        },
        library_name_labels => {
            $self->get_library_name_labels('tumor', $self->tumor_sample, [$self->build]),
            $self->get_library_name_labels('normal', $self->normal_sample, [$self->build]),
        },
    };
}

sub get_sample_name_labels {
    my ($self, $category) = @_;

    my $accessor = sprintf('%s_sample', $category);
    return (
        $self->model->$accessor->name  => sprintf('%s(%s)', ucfirst($category), $self->model->$accessor->name),
    );
}

sub is_valid {
    my $self = shift;

    if (my @problems = $self->__errors__) {
        $self->error_message('Germline is invalid!');
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

sub is_single_bam {
    my $self = shift;
    my $model = shift;
    return (!defined($model->normal_sample));
}

sub build {
    my $self = shift;
    return $self->model->last_succeeded_build;
}

1;

