package Genome::TestObjGenerator::Model;
use Genome::TestObjGenerator::Base;
@ISA = (Genome::TestObjGenerator::Base);

use strict;
use warnings;
use Genome;
use Genome::TestObjGenerator::Util;

our @required_params = qw(name processing_profile_id);

sub create_name {
    return Genome::TestObjGenerator::Util::generate_name("test_model");
}

1;
