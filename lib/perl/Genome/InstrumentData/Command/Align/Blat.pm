package Genome::InstrumentData::Command::Align::Blat;

#REVIEW fdu
#limited use, removable, see REVIEW in base class Align.pm

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Align::Blat {
    is => ['Genome::InstrumentData::Command::Align'],
    has_constant => [
        aligner_name                    => { value => 'blat' },
    ],
    has_param => [
        version                 => { default_value => '34'},
    ],
    doc => 'align instrument data using Blat'
};

sub help_synopsis {
return <<EOS
FIXME
EOS
}

sub help_detail {
return <<EOS
Launch the Blat aligner in a standard way and produce results ready for the genome modeling pipeline.

EOS
}


1;

