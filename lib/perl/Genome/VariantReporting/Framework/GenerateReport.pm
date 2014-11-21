package Genome::VariantReporting::Framework::GenerateReport;

use strict;
use warnings;
use Genome;
use Try::Tiny qw(try catch);
use Digest::MD5 qw(md5_hex);

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
        variant_type => {
            is => 'Text',
            valid_values => ['snvs', 'indels'],
        },
        provider_json => {
            is => 'Text',
        },
        user => {
            is => 'Genome::Process',
            is_optional => 1,
            id_by => 'process_id',
        },
        process_id => {
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
        plan_json_lookup => md5_hex($self->plan_json),
        variant_type => $self->variant_type,
        provider_json => $self->provider_json,
        provider_json_lookup => md5_hex($self->provider_json),
        test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME},
    );
}


1;
