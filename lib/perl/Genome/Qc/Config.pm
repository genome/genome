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
        PicardCollectMultipleMetrics => {
            params => {
                param1 => 'a',
                param2 => 'b',
            },
            dependencies => [],
        },
    };
}

1;

