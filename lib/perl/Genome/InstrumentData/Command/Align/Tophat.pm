package Genome::InstrumentData::Command::Align::Tophat;

#REVIEW fdu
#limited use, removable, see REVIEW in base class Align.pm

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Align::Tophat {
    is => ['Genome::InstrumentData::Command::Align'],
    has_constant => [
        aligner_name => { value => 'tophat' },
    ],
    doc => "align instrument data using Tophat which in turn uses Bowtie(http://tophat.cbcb.umd.edu/)",
};

sub help_synopsis {
return <<EOS
genome instrument-data align tophat -r NCBI-human-build36 -i 2761701954

genome instrument-data align tophat -r NCBI-human-build36 -i 2761701954 -v 2.03.01 

genome instrument-data align tophat --reference-name NCBI-human-build36 --instrument-data-id 2761701954 --version 2.03.01 

genome instrument-data align tophat -i 2761701954 -v 2.03.01 
EOS
}

sub help_detail {
return <<EOS
Launch the Tophat aligner in a standard way and produce results ready for the genome modeling pipeline.

EOS
}


1;

