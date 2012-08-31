package Genome::InstrumentData::Command::Align::Bwa;

#REVIEW fdu
#limited use, removable, see REVIEW in base class Align.pm

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Align::Bwa {
    is => ['Genome::InstrumentData::Command::Align'],
    has_constant => [
        aligner_name                    => { value => 'bwa' },
    ],
    has_param => [
        version                 => { default_value => '0.5.5'},
    ],
    doc => 'align instrument data using BWA (see http://maq.sourceforge.net)',
};

sub help_synopsis {
return <<EOS
FIXME
EOS
}

sub help_detail {
return <<EOS
Launch the bwa aligner in a standard way and produce results ready for the genome modeling pipeline.

See http://maq.sourceforge.net.
EOS
}


1;

