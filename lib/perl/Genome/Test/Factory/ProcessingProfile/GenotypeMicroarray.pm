package Genome::Test::Factory::ProcessingProfile::GenotypeMicroarray;
use Genome::Test::Factory::ProcessingProfile;
@ISA = (Genome::Test::Factory::ProcessingProfile);

use strict;
use warnings;

sub create_id {
    return Genome::Test::Factory::Util::generate_name("pp_id");
}

1;
