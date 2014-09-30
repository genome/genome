package Genome::Test::Factory::ProcessingProfile::SomaticVariation;
use Genome::Test::Factory::ProcessingProfile;
@ISA = (Genome::Test::Factory::ProcessingProfile);

use strict;
use warnings;

our @required_params = qw(tiering_version bam_readcount_version);

sub create_tiering_version {
    return 1;
}

sub create_bam_readcount_version {
    return 0.6;
}

1;
