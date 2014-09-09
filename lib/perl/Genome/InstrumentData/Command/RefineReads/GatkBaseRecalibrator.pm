package Genome::InstrumentData::Command::RefineReads::GatkBaseRecalibrator;

use strict;
use warnings;

use Genome;

use Genome::InstrumentData::Gatk::BaseRecalibratorBamResult;

class Genome::InstrumentData::Command::RefineReads::GatkBaseRecalibrator {
    is => 'Genome::InstrumentData::Command::RefineReads::GatkBase',
};
Genome::InstrumentData::Command::RefineReads::GatkBaseRecalibrator->__meta__->property_meta_for_name('known_sites')->is_optional(undef);

sub result_names {
    return ( 'base recalibrator bam' );
}

1;

