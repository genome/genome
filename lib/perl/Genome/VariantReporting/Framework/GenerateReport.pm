package Genome::VariantReporting::Framework::GenerateReport;

use strict;
use warnings;
use Genome;
use Try::Tiny qw(try catch);

use Genome::VariantReporting::Framework::FileLookup qw(
    calculate_lookup
);

class Genome::VariantReporting::Framework::GenerateReport {
    is => 'Genome::Command::DelegatesToResult',
    has_input => [
        input_vcf => {
            is => 'File',
        },
        plan_json => {
            is => 'Text',
        },
        output_directory => {
            is => 'Path',
            is_output => 1,
        },
        variant_type => {
            is => 'Text',
            valid_values => ['snvs', 'indels'],
        },
        provider_json => {
            is => 'Text',
        },
    ],
    has_param => [
        lsf_resource => {
            value => q{-R 'select[mem>16000] rusage[mem=16000]' -M 16000000},
        },
    ],
};

sub result_class {
    return 'Genome::VariantReporting::Framework::ReportResult';
}

sub input_hash {
    my $self = shift;
    return (
        input_vcf => $self->input_vcf,
        input_vcf_lookup => calculate_lookup($self->input_vcf),
        plan_json => $self->plan_json,
        variant_type => $self->variant_type,
        provider_json => $self->provider_json,
        test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME},
    );
}

sub post_get_or_create {
    my $self = shift;
    $self->symlink_results;
    return 1;
}

sub symlink_results {
    my $self = shift;

    try {
        Genome::Sys->create_directory($self->output_directory);
        Genome::Sys->symlink_directory(
            $self->output_result->output_dir,
            $self->output_directory,
        );
    } catch {
        $self->error_message("Could not symlink to output_directory because %s", $_);
        $self->output_directory($self->output_result->output_dir);
    };
    $self->status_message("Outputs located at %s", $self->output_directory);
}


1;
