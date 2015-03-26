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
        picard_collect_multiple_metrics => {
            class => 'Genome::Qc::Tool::PicardCollectMultipleMetrics',
            params => {
                param1 => 'a',
                param2 => 'b',
            },
            dependency => {name => "bam_file", fd => "STDOUT"},
        },
        samtools_view => {
            class => 'Genome::Qc::Tool::BamFile',
            params => {
                input_file => '-',
            },
            in_file => "bam_file",
        },
    };
}

1;

