package Genome::InstrumentData::Command::Align::Brat;

#REVIEW fdu
#limited use, removable, see REVIEW in base class Align.pm

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Align::Brat {
    is => ['Genome::InstrumentData::Command::Align'],
    has_constant => [
        aligner_name                    => { value => 'brat' },
    ],
    has_param => [
        version                 => { default_value => '1.2.1-mod'},
    ],
    doc => 'align instrument data using BRAT (see http://compbio.cs.ucr.edu/brat/)',
};

sub help_synopsis {
return <<EOS
FIXME
EOS
}

sub help_detail {
return <<EOS
Launch the brat aligner in a standard way and produce results ready for the genome modeling pipeline.

See http://compbio.cs.ucr.edu/brat/.
EOS
}


1;

