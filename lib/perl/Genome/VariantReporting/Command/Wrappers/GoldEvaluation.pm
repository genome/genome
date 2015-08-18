package Genome::VariantReporting::Command::Wrappers::GoldEvaluation;

use strict;
use warnings;

use Genome;

use File::Basename qw(dirname);
use File::Spec;

class Genome::VariantReporting::Command::Wrappers::GoldEvaluation {
    is => ['Command::V2', 'Genome::VariantReporting::Command::Wrappers::Utils'],
    has_input => {
        model => {
            is => 'Genome::Model::SomaticValidation',
        },
        gold_sample_name => {
            is => 'Text',
            doc => 'The name of the sample in the VCF used to retrieve genotypes.',
        },
        snvs_input_vcf => {
            is => 'Path',
        },
        indels_input_vcf => {
            is => 'Path',
        },
    },

    has => [
        normal_sample => {
          is => 'Genome::Sample',
          via => 'model',
          to => 'normal_sample',
        },

        normal_sample_name => {
          is => 'Text',
          via => 'normal_sample',
          to => 'name',
        },

        discovery_sample => {
          is => 'Genome::Sample',
          via => 'model',
          to => 'tumor_sample',
        },

        discovery_sample_name => {
          is => 'Text',
          via => 'discovery_sample',
          to => 'name',
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
            common_translations => $self->get_germline_translations(),
            discovery => $self->build,
            plan_file_basename => 'gold_germline_report_TYPE.yaml',
            gold_sample_name => $self->gold_sample_name,
            label => 'gold_germline',
        );
    } else {
        #Somatic
        $model_pair = Genome::VariantReporting::Command::Wrappers::ModelPair->create(
            common_translations => $self->get_somatic_translations(),
            discovery => $self->build,
            plan_file_basename => 'gold_somatic_report_TYPE.yaml',
            gold_sample_name => $self->gold_sample_name,
            label => 'gold_somatic',
        );
    }
    for my $variant_type (qw(snvs indels)) {
        my $vcf_method = $variant_type .'_input_vcf';
        my %params = (
            input_vcf => $self->$vcf_method,
            variant_type => $variant_type,
            plan_file => $model_pair->plan_file($variant_type),
            translations_file => $model_pair->translations_file,
        );
        Genome::VariantReporting::Command::CreateReport->execute(%params);
    }
    return 1;
};

sub build {
    my $self = shift;
    return $self->model->last_succeeded_build;
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

sub get_somatic_translations {
    my $self = shift;

    return {
        sample_name_labels => {
            $self->discovery_sample_name =>
                sprintf('Discovery(%s)', $self->discovery_sample_name),
            $self->normal_sample_name =>
                sprintf('Normal(%s)', $self->normal_sample_name),
            $self->gold_sample_name =>
                sprintf('Gold(%s)', $self->gold_sample_name),
        },
        library_name_labels => {
            $self->get_library_name_labels('discovery'),
            $self->get_library_name_labels('normal'),
        },
    };
}

sub get_germline_translations {
    my $self = shift;

    return {
        sample_name_labels => {
            $self->discovery_sample_name =>
                sprintf('Discovery(%s)', $self->discovery_sample_name),
            $self->gold_sample_name =>
                sprintf('Gold(%s)', $self->gold_sample_name),
        },
        library_name_labels => {
            $self->get_library_name_labels('discovery'),
        },
    };
}

1;
