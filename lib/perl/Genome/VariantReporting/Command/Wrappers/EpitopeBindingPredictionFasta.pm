package Genome::VariantReporting::Command::Wrappers::EpitopeBindingPredictionFasta;

use strict;
use warnings;
use Genome;
use File::Basename qw(dirname);
use YAML;

class Genome::VariantReporting::Command::Wrappers::EpitopeBindingPredictionFasta {
    is => 'Command::V2',
    has_input => [
        build => {
            is => "Genome::Model::Build::SomaticVariation",
        },
        variant_type => {
            is => 'Text',
            valid_values => ['indels', 'snvs'],
        },
    ],
};

sub execute {
    my $self = shift;

    $self->_validate_inputs;

    my $report = $self->_create_report;
    $report->execute;
}

sub _validate_inputs {
    my $self = shift;

    unless ($self->has_valid_variant_type_for_build) {
        $self->fatal_message('Input vcf of type (%s) does not exist for build (%s)', $self->variant_type, $self->build->id);
    }
}

sub has_valid_variant_type_for_build {
    my $self = shift;

    my $input_vcf = $self->input_vcf($self->variant_type);
    return (-s $input_vcf);
}

sub _create_report {
    my $self = shift;

    my $variant_type = $self->variant_type;
    my %params = (
        input_vcf => $self->input_vcf($variant_type),
        variant_type => $variant_type,
        plan_file => $self->plan_file($variant_type),
        translations_file => $self->translations_file,
    );

    return Genome::VariantReporting::Command::CreateReport->create(%params);
}

sub report_workflow {
    my $self = shift;

    $self->_create_report->dag;
}

sub input_vcf {
    my ($self, $variant_type) = @_;
    my $accessor = "get_detailed_${variant_type}_vcf";
    return $self->build->$accessor;
}

sub plan_file {
    my ($self, $variant_type) = @_;
    my $base_dir = Genome->base_dir;
    return File::Spec->join($base_dir, 'VariantReporting', 'plan_files', "epitope_prediction_$variant_type.yaml");
}

sub translations_file {
    my $self = shift;
    my $translations_file = Genome::Sys->create_temp_file_path;
    my %translations;
    $translations{reference_fasta} = $self->build->reference_sequence_build->full_consensus_path("fa");
    YAML::DumpFile(File::Spec->join($translations_file), \%translations);
    return $translations_file;
}

1;
