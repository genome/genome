package Genome::Qc::Config;

use strict;
use warnings;
use Genome;

class Genome::Qc::Config {
    has => [
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
        bam_file => {
            class => 'Genome::Qc::Tool::BamFile',
        },
    };
}

1;

