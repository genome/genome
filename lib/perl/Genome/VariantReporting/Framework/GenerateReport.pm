package Genome::VariantReporting::Framework::GenerateReport;

use strict;
use warnings;
use Genome;
use Data::Dump qw(pp);

use Genome::VariantReporting::Framework::FileLookup qw(
    calculate_lookup
);

class Genome::VariantReporting::Framework::GenerateReport {
    is => 'Command::V2',
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
    has_optional_output => [
        software_result => {
            is => 'Genome::SoftwareResult',
        },
    ],
};

sub result_class {
    return 'Genome::VariantReporting::Framework::ReportResult';
}

sub shortcut {
    my $self = shift;

    $self->debug_message("Attempting to get a %s with arugments %s",
        $self->result_class, pp($self->input_hash));
    my $result = $self->result_class->get_with_lock($self->input_hash);
    if ($result) {
        $self->debug_message("Found existing result (%s)", $result->id);
        $self->software_result($result);
        $self->symlink_results;
        return 1;
    } else {
        $self->debug_message("Found no existing result.");
        return 0;
    }
}

sub execute {
    my $self = shift;

    $self->debug_message("Attempting to get or create a %s with arugments %s",
        $self->result_class, pp({$self->input_hash}));
    my $result = $self->result_class->get_or_create($self->input_hash);
    $self->debug_message("Got or created result (%s)", $result->id);
    $self->debug_message("Calculated query for result:\n".pp({$result->calculate_query()}));
    $self->software_result($result);
    $self->symlink_results;
    return 1;
}

sub symlink_results {
    my $self = shift;

    local $@;
    eval {
        Genome::Sys->create_directory($self->output_directory);
        Genome::Sys->symlink_directory(
            $self->software_result->output_dir,
            $self->output_directory,
        );
    };
    my $error = $@;
    if ($error) {
        $self->error_message("Could not symlink to output_directory because %s", $error);
        $self->output_directory($self->software_result->output_dir);
    }
    $self->status_message("Outputs located at %s", $self->output_directory);
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


1;
