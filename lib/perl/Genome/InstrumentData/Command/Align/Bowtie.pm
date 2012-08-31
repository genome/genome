package Genome::InstrumentData::Command::Align::Bowtie;

#REVIEW fdu
#limited use, removable, see REVIEW in base class Align.pm

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Align::Bowtie {
    is => ['Genome::InstrumentData::Command::Align'],
    has_constant => [
        aligner_name => { value => 'Bowtie' },
    ],
    has => [
        version => {is=>'String', default_value=>Genome::Model::Tools::Bowtie->default_bowtie_version}
    ],
    doc => "align instrument data using Bowtie's novoalign tool (see http://bowtie-bio.sourceforge.net)",
};

sub help_synopsis {
# TODO: Make these actually reflect real examples.
return <<EOS
genome instrument-data align bowtie -r NCBI-human-build36 -i 2761701954

genome instrument-data align bowtie -r NCBI-human-build36 -i 2761701954 -v 2.03.01 

genome instrument-data align bowtie --reference-name NCBI-human-build36 --instrument-data-id 2761701954 --version 2.03.01 

genome instrument-data align bowtie -i 2761701954 -v 2.03.01 
EOS
}

sub help_detail {
return <<EOS
Launch the Bowtie aligner in a standard way and produce results ready for the genome modeling pipeline.

See http://bowtie-bio.sourceforge.net.
EOS
}


1;

