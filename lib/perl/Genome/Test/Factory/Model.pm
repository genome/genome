package Genome::Test::Factory::Model;
use Genome::Test::Factory::Base;
@ISA = (Genome::Test::Factory::Base);

use strict;
use warnings;
use Genome;
use Genome::Test::Factory::Util;

our @required_params = qw(name processing_profile_id);

sub create_name {
    return Genome::Test::Factory::Util::generate_name("test_model");
}

1;
