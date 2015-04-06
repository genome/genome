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

