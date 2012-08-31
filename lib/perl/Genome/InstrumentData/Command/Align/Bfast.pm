package Genome::InstrumentData::Command::Align::Bfast;

#REVIEW fdu
#limited use, removable, see REVIEW in base class Align.pm

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Align::Bfast {
    is => ['Genome::InstrumentData::Command::Align'],
    has_constant => [
        aligner_name                    => { value => 'bfast' },
    ],
    has_param => [
        version                 => { default_value => '0.6.4d'},
    ],
    doc => 'align instrument data using Bfast'
};

sub help_synopsis {
return <<EOS
FIXME
EOS
}

sub help_detail {
return <<EOS
Launch the Bfast aligner in a standard way and produce results ready for the genome modeling pipeline.

EOS
}


1;

