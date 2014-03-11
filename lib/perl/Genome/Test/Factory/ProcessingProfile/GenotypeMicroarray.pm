package Genome::Test::Factory::ProcessingProfile::GenotypeMicroarray;
use Genome::Test::Factory::ProcessingProfile;
@ISA = (Genome::Test::Factory::ProcessingProfile);

use strict;
use warnings;

our @required_params = qw(instrument_type id);

sub create_instrument_type {
    return "unknown";
}

sub create_id {
    return Genome::Test::Factory::Util::generate_name("pp_id");
}
1;
