package Genome::InstrumentData::Command::Align::Shrimp2;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Align::Shrimp2 {
    is => ['Genome::InstrumentData::Command::Align'],
    has_constant => [
        aligner_name                    => { value => 'shrimp2' },
    ],
    has_param => [
        version                 => { default_value => '2.0.1'},
    ],
    doc => 'align instrument data using SHRiMP2'
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

