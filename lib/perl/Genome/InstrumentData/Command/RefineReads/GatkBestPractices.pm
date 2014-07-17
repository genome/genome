package Genome::InstrumentData::Command::RefineReads::GatkBestPractices;

use strict;
use warnings;

use Genome;

use Genome::InstrumentData::Gatk::IndelRealignerResult;
use Genome::InstrumentData::Gatk::BaseRecalibratorBamResult;
require Genome::Utility::Text;

class Genome::InstrumentData::Command::RefineReads::GatkBestPractices {
    is => 'Genome::InstrumentData::Command::RefineReads::GatkBase',
    has_optional_calculated => {
        indel_realigner_result => { calculate => q( return $self->results->{indel_realigner_result}; ), },
        base_recalibrator_bam_result => { calculate => q( return $self->results->{base_recalibrator_bam_result}; ), },
    },
};

sub result_names {
    return ( 'indel realigner', 'base recalibrator bam' );
}

1;

