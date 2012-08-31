package Genome::InstrumentData::Command::Align::Smalt;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Align::Smalt {
    is => ['Genome::InstrumentData::Command::Align'],
    has_constant => [
        aligner_name => { value => 'smalt' },
    ],
    has_param => [
        version => { default_value => '0.5.5'},
    ],
    doc => 'align instrument data using smalt'
};

sub help_synopsis {
return <<EOS
FIXME
EOS
}

sub help_detail {
return <<EOS
Launch the aligner in a standard way and produce results ready for the genome modeling pipeline.

EOS
}


1;

