package Genome::InstrumentData::Command::RefineReads::GatkIndelRealigner;

use strict;
use warnings;

use Genome;

use Genome::InstrumentData::Gatk::IndelRealignerResult;

class Genome::InstrumentData::Command::RefineReads::GatkIndelRealigner {
    is => 'Genome::InstrumentData::Command::RefineReads::GatkBase',
};

sub result_names {
    return ( 'indel realigner' );
}

1;

