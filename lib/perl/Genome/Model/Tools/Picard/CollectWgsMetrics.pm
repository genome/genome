package Genome::Model::Tools::Picard::CollectWgsMetrics;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Picard::CollectWgsMetrics {
    is => 'Genome::Model::Tools::Picard::Base',
    has_input => [
        input_file   => {
            is  => 'String',
            doc => 'The SAM/BAM files to run on.  File type is determined by suffix.',
            picard_param_name => 'INPUT',
        },
        output_file  => {
            is  => 'String',
            doc => 'The output metrics file',
            picard_param_name => 'OUTPUT',
        },
        reference_sequence => {
            is => 'String',
            doc => 'The reference sequence fasta aligned to.',
            picard_param_name => 'REFERENCE_SEQUENCE',
        },
        minimum_mapping_quality => {
            is => 'Integer',
            is_optional => 1,
            doc => "Minimum mapping quality for a read to contribute coverage. Default value: 20. This option can be set to 'null' to clear the default value.",
            picard_param_name => 'MINIMUM_MAPPING_QUALITY',
        },
        minimum_base_quality => {
            is => 'Integer',
            is_optional => 1,
            doc => "Minimum base quality for a base to contribute coverage. Default value: 20. This option can be set to 'null' to clear the default value.",
            picard_param_name => 'MINIMUM_BASE_QUALITY',
        },
        coverage_cap => {
            is => 'Integer',
            is_optional => 1,
            doc => "Treat bases with coverage exceeding this value as if they had coverage at this value. Default value: 250. This option can be set to 'null' to clear the default value.",
            picard_param_name => 'COVERAGE_CAP',
        },
        stop_after => {
            is => 'Integer',
            is_optional => 1,
            doc => "For debugging purposes, stop after processing this many genomic bases. Default value: -1. This option can be set to 'null' to clear the default value.",
            picard_param_name => 'STOP_AFTER',
        },
        include_bq_histogram => {
            is => 'Boolean',
            is_optional => 1,
            doc => "Determines whether to include the base quality histogram in the metrics file.",
            picard_param_name => 'INCLUDE_BQ_HISTORAM',
        },
    ],
};

sub help_brief {
    'Computes a number of metrics that are useful for evaluating coverage and performance of whole genome sequencing experiments.';
}

sub help_detail {
    return <<EOS
    For Picard documentation of this command see:
    http://broadinstitute.github.io/picard/command-line-overview.html#CollectWgsMetrics
EOS
}

sub _jar_name {
    return 'CollectWgsMetrics.jar';
}

sub _java_class {
    return qw(picard analysis CollectWgsMetrics);
}

sub _validate_params {
    my $self = shift;
    $self->enforce_minimum_version('1.114');
}

1;
