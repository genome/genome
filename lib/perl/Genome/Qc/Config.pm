package Genome::Qc::Config;

use strict;
use warnings;
use Genome;

class Genome::Qc::Config {
    is => 'UR::Value',
    id_by => [
        name => {
            is => 'String',
        },
    ],
};

sub get_commands_for_alignment_result {
    return {
        picard_calculate_hs_metrics => {
            class => 'Genome::Qc::Tool::Picard::CalculateHsMetrics',
            params => {
                input_file => '/dev/stdin',
                bait_intervals => 'bait_intervals', #region_of_interest_set
                target_intervals => 'target_intervals', #target_region_set
                use_version => 1.123,
            },
            in_file => "bam_file",
        },
        picard_collect_wgs_metrics => {
            class => 'Genome::Qc::Tool::Picard::CollectWgsMetrics',
            params => {
                input_file => '/dev/stdin',
                reference_sequence => 'reference_sequence',
                use_version => 1.123,
            },
            in_file => "bam_file",
        },
        picard_collect_gc_bias_metrics => {
            class => 'Genome::Qc::Tool::Picard::CollectGcBiasMetrics',
            params => {
                input_file => '/dev/stdin',
                refseq_file => 'reference_sequence',
                assume_sorted => 1,
                use_version => 1.123,
                output_file=> '/dev/null',
            },
            in_file => "bam_file",
        },
        picard_mark_duplicates => {
            class => 'Genome::Qc::Tool::Picard::MarkDuplicates',
            params => {
                output_file => '/dev/null',
                input_file => '/dev/stdin',
                use_version => 1.123,
            },
            in_file => "bam_file",
        },
        picard_collect_multiple_metrics => {
            class => 'Genome::Qc::Tool::Picard::CollectMultipleMetrics',
            params => {
                input_file => '/dev/stdin',
                use_version => 1.123,
            },
            in_file => "bam_file",
        },
    };
}

1;

