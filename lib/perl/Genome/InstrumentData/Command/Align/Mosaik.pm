package Genome::InstrumentData::Command::Align::Mosaik;

#REVIEW fdu
#limited use, removable, see REVIEW in base class Align.pm

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Align::Mosaik {
    is => ['Genome::InstrumentData::Command::Align'],
    has_constant => [
        aligner_name                    => { value => 'mosaik' },
    ],
    has_param => [
        version                 => { default_value => '1.0.1388'},
    ],
    doc => 'align instrument data using Mosaik'
};

sub help_synopsis {
return <<EOS
FIXME
EOS
}

sub help_detail {
return <<EOS
Launch the Mosaik aligner in a standard way and produce results ready for the genome modeling pipeline.

EOS
}


1;

