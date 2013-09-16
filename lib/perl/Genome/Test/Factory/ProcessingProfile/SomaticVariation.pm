package Genome::Test::Factory::ProcessingProfile::SomaticVariation;
use Genome::Test::Factory::ProcessingProfile;
@ISA = (Genome::Test::Factory::ProcessingProfile);

use strict;
use warnings;

our @required_params = qw(tiering_version);

sub create_tiering_version {
    return 1;
}

1;
