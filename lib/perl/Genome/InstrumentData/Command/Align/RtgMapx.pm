package Genome::InstrumentData::Command::Align::RtgMapx;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Align::RtgMapx {
    is => ['Genome::InstrumentData::Command::Align'],
    has_constant => [
        aligner_name                    => { value => 'rtg mapx' },
    ],
    has_param => [
        version                 => { default_value => 'v2.0.1-build-28762'},
    ],
    doc => 'align instrument data using RTG (see ~rtg/rtg/README.txt)',
};

sub help_synopsis {
return <<EOS
FIXME
EOS
}

sub help_detail {
return <<EOS
Launch the rtp map aligner in a standard way and produce results ready for the genome modeling pipeline.

See ~rtg/rtg/README.txt
EOS
}


1;

