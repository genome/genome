package Genome::InstrumentData::Command::Align::Crossmatch;

#REVIEW fdu
#limited use, removable, see REVIEW in base class Align.pm

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Align::Crossmatch {
    is => ['Genome::InstrumentData::Command::Align'],
    has_constant => [
        aligner_name                    => { value => 'crossmatch' },
    ],
    has_param => [
        version                 => { default_value => '1.080721'},
    ],
    doc => 'align instrument data using Crossmatch'
};

sub help_synopsis {
return <<EOS
FIXME
EOS
}

sub help_detail {
return <<EOS
Launch the Crossmatch aligner in a standard way and produce results ready for the genome modeling pipeline.

EOS
}


1;

