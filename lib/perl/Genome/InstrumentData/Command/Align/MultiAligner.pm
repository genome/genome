package Genome::InstrumentData::Command::Align::MultiAligner;

#REVIEW fdu
#limited use, removable, see REVIEW in base class Align.pm

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Align::MultiAligner {
    is => ['Genome::InstrumentData::Command::Align'],
    has_constant => [
        aligner_name                    => { value => 'multi aligner' },
    ],
    has_param => [
        version                 => { default_value => '0'},
    ],
    doc => 'align instrument data using multiple aligners',
};

sub help_synopsis {
return <<EOS
FIXME
EOS
}

sub help_detail {
return <<EOS
Launch a combination of aligners in a standard way and produce results ready for the genome modeling pipeline.
EOS
}


1;

