package Genome::InstrumentData::Command::RefineReads::GatkIndelRealigner;

use strict;
use warnings;

use Genome;

use Genome::InstrumentData::Gatk::IndelRealignerResult;

class Genome::InstrumentData::Command::RefineReads::GatkIndelRealigner {
    is => 'Genome::InstrumentData::Command::RefineReads::GatkBase',
    has_optional_calculated => {
        indel_realigner_result => { calculate => q( return $self->results->{indel_realigner_result}; ), },
    },
};

sub result_names {
    return ( 'indel realigner' );
}

1;

