package Genome::InstrumentData::Command::RefineReads::GatkBaseRecalibrator;

use strict;
use warnings;

use Genome;

use Genome::InstrumentData::Gatk::BaseRecalibratorBamResult;

class Genome::InstrumentData::Command::RefineReads::GatkBaseRecalibrator {
    is => 'Genome::InstrumentData::Command::RefineReads::GatkBase',
    has_optional_calculated => {
        base_recalibrator_bam_result => { calculate => q( return $self->results->{base_recalibrator_bam_result}; ), },
    },
};

sub result_names {
    return ( 'base recalibrator bam' );
}

1;

