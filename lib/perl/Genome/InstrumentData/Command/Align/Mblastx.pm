package Genome::InstrumentData::Command::Align::Mblastx;

use strict;
use warnings;
use Genome;

class Genome::InstrumentData::Command::Align::Mblastx {
    is => ['Genome::InstrumentData::Command::Align'],
    has_constant => [
        aligner_name  => { 
            value => 'mblastx', 
        },
    ],
    has_param => [
        version => { 
            default_value => '09242010', 
        },
        mhashgen_format => { 
            is => 'Text', 
            valid_values=>[
                "K",
                "N",
                "A",
            ],
        },
    ],
    doc => 'align instrument data using MCW mblastx',
};

sub help_synopsis {
return <<EOS
FIXME
EOS
}

sub help_detail {
return <<EOS
Launch the MCW mblastx aligner in a standard way and produce results ready for the genome modeling pipeline.

EOS
}

1;

