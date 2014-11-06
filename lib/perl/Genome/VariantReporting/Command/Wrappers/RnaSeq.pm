package Genome::VariantReporting::Command::Wrappers::RnaSeq;

use strict;
use warnings;

use Genome;

use File::Basename qw(dirname);
use File::Spec;
use Genome::Utility::Text;

class Genome::VariantReporting::Command::Wrappers::RnaSeq {
    is => 'Command::V2',
    has => {
        somatic_build => {
            is => 'Genome::Model::Build',
        },
        tumor_build => {
            is => 'Genome::Model::Build::RnaSeq',
            doc => 'The build that contains the fpkm file',
        },
        base_output_dir => { is => 'Text', },
    },
    has_calculated_optional => {
        output_dir => {
            calculate_from => [qw/ base_output_dir somatic_build/],
            calculate => q| my $model_name  = $somatic_build->model->name;
                my $safe_model_name = Genome::Utility::Text::sanitize_string_for_filesystem($model_name);
                return File::Spec->join($base_output_dir, $safe_model_name); |,
        },
        translations_file => {
            calculate_from => [qw/ output_dir /],
            calculate => q( File::Spec->join($output_dir, "resource.yaml") ),
        },
    },
};

sub reports_directory {
    my ($self, $variant_type) = @_;
    return File::Spec->join($self->output_dir, "reports_$variant_type");
};

sub logs_directory {
    my ($self, $variant_type) = @_;
    return  File::Spec->join($self->output_dir, "logs_$variant_type");
};

sub execute {
    my $self = shift;
    Genome::Sys->create_directory($self->output_dir);
    Genome::Sys->create_directory($self->reports_directory("snvs"));
    Genome::Sys->create_directory($self->logs_directory("snvs"));
    $self->generate_translations_file;
    $self->run_reports;
    return 1;
};

sub is_valid {
    my $self = shift;

    if (my @problems = $self->__errors__) {
        $self->error_message('RnaSeq is invalid!');
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

sub generate_translations_file {
    my $self = shift;

    return if not $self->is_valid;
    my %translations;

    my @aligned_bam_results;
    if ($self->somatic_build->isa('Genome::Model::Build::SomaticValidation')) {
        push @aligned_bam_results, $self->somatic_build->merged_alignment_result->id;
        push @aligned_bam_results, $self->somatic_build->control_merged_alignment_result->id;
    }
    elsif ($self->somatic_build->isa('Genome::Model::Build::SomaticVariation')) {
        push @aligned_bam_results, $self->somatic_build->tumor_build->merged_alignment_result->id;
        push @aligned_bam_results, $self->somatic_build->normal_build->merged_alignment_result->id;
    }
    else {
        die $self->error_message("somatic_build is of unhandled type: (%s). Needs to be either 'Genome::Model::Build::SomaticValidation' or 'Genome::Model::Build::SomaticVariation'", $self->somatic_build->class);
    }
    $translations{aligned_bam_result_id} = \@aligned_bam_results;

    $translations{reference_fasta} = $self->somatic_build->reference_sequence_build->full_consensus_path("fa");

    $translations{feature_list_ids} = {};

    $translations{fpkm_file} = File::Spec->join($self->tumor_build->data_directory, 'expression', 'genes.fpkm_tracking');

    if ($self->somatic_build->isa('Genome::Model::Build::SomaticValidation')) {
        $translations{tumor} = $self->somatic_build->tumor_sample->name;
        $translations{normal} = $self->somatic_build->normal_sample->name;
    }
    elsif ($self->somatic_build->isa('Genome::Model::Build::SomaticVariation')) {
        $translations{tumor} = $self->somatic_build->tumor_build->subject->name;
        $translations{normal} = $self->somatic_build->normal_build->subject->name;
    }

    YAML::DumpFile($self->translations_file, \%translations);

    return 1;
}

sub run_reports {
    my $self = shift;

    my $variant_type = 'snvs';
    Genome::VariantReporting::Command::CreateReport->execute(
        input_vcf => $self->input_vcf($variant_type),
        variant_type => $variant_type,
        output_directory => $self->reports_directory($variant_type),
        plan_file => $self->plan_file($variant_type),
        translations_file => $self->translations_file,
        log_directory => $self->logs_directory($variant_type),
    );

}

sub plan_file {
    my ($self, $type) = @_;
    return File::Spec->join($self->_plan_search_dir, "rnaseq_variants_$type.yaml");
}

sub input_vcf {
    my ($self, $variant_type) = @_;
    return $self->somatic_build->get_detailed_vcf_result($variant_type)->get_vcf($variant_type);
}

sub _plan_search_dir {
    my $variant_reporting_base_dir = dirname(dirname(dirname(__FILE__)));
    return File::Spec->join($variant_reporting_base_dir, 'plan_files');
}

1;

