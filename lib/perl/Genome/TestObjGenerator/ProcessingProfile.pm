package Genome::TestObjGenerator::ProcessingProfile;
use Genome::TestObjGenerator::Base;
@ISA = (Genome::TestObjGenerator::Base);

use strict;
use warnings;
use Genome;
use Genome::TestObjGenerator::Util;

my @required_params = ("name");

sub generate_obj {
    my $class = shift;
    (my $pp_class = $class) =~ s/::TestObjGenerator::/::/;
    return $pp_class->create(@_);
}

sub get_required_params {
    return \@required_params;
}

sub create_name {
    return Genome::TestObjGenerator::Util::generate_name("test_processing_profile");
}

1;
