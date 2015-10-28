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
        output_directory => {
            is => 'Path',
        },
    ],
    has_calculated_optional => [
        translations_file => {
            calculate_from => [qw(output_directory)],
            calculate => q/File::Spec->join($output_directory, "resources.yaml")/,
        },
    ],
};

sub execute {
    my $self = shift;
    $self->generate_translations_file;
    $self->run_reports;
}


sub generate_translations_file {
    my $self = shift;
    my %translations;
    $translations{reference_fasta} = $self->build->reference_sequence_build->full_consensus_path("fa");
    YAML::DumpFile(File::Spec->join($self->translations_file), \%translations);
}

sub run_reports {
    my $self = shift;
    for my $variant_type (qw(snvs indels)) {
        my $input_vcf = $self->input_vcf($variant_type);
        next unless -s $input_vcf;

        my %params = (
            input_vcf => $input_vcf,
            variant_type => $variant_type,
            plan_file => $self->plan_file($variant_type),
            translations_file => $self->translations_file,
        );
        Genome::VariantReporting::Command::CreateReport->execute(%params);
    }
}

sub plan_file {
    my ($self, $variant_type) = @_;
    my $base_dir = dirname(dirname(dirname(__FILE__)));
    return File::Spec->join($base_dir, "plan_files", "epitope_prediction_$variant_type.yaml");
}

sub input_vcf {
    my ($self, $variant_type) = @_;
    my $accessor = "get_detailed_$variant_type"."_vcf";
    return $self->build->$accessor;
}

sub report_directory {
    my ($self, $variant_type) = @_;
    return File::Spec->join($self->output_directory, "report_$variant_type");
}

1;
