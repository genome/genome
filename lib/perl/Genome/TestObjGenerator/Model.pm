package Genome::TestObjGenerator::Model;
use Genome::TestObjGenerator::Base;
@ISA = (Genome::TestObjGenerator::Base);

use strict;
use warnings;
use Genome;
use Genome::TestObjGenerator::Util;

my @required_params = ("name", "processing_profile_id");

sub get_required_params {
    return \@required_params;
}

sub create_name {
    return Genome::TestObjGenerator::Util::generate_name("test_model");
}

1;

