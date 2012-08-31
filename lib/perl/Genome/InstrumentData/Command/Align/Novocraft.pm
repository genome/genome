package Genome::InstrumentData::Command::Align::Novocraft;

#REVIEW fdu
#limited use, removable, see REVIEW in base class Align.pm

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Align::Novocraft {
    is => ['Genome::InstrumentData::Command::Align'],
    has_constant => [
        aligner_name => { value => 'novocraft' },
    ],
    has => [
        version => {is=>'String', default_value=>Genome::Model::Tools::Novocraft->default_novocraft_version}
    ],
    doc => "align instrument data using novocraft's novoalign tool (see http://novocraft.com)",
};

sub help_synopsis {
return <<EOS
genome instrument-data align novocraft -r NCBI-human-build36 -i 2761701954

genome instrument-data align novocraft -r NCBI-human-build36 -i 2761701954 -v 2.03.01 

genome instrument-data align novocraft --reference-name NCBI-human-build36 --instrument-data-id 2761701954 --version 2.03.01 

genome instrument-data align novocraft -i 2761701954 -v 2.03.01 
EOS
}

sub help_detail {
return <<EOS
Launch the novocraft novalign aligner in a standard way and produce results ready for the genome modeling pipeline.

See http://novocraft.com.
EOS
}


1;

