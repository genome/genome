package Genome::TestObjGenerator::ProcessingProfile::SomaticVariation;
use Genome::TestObjGenerator::ProcessingProfile;
@ISA = (Genome::TestObjGenerator::ProcessingProfile);

use strict;
use warnings;
use Genome;

my @required_params = ("tiering_version");

sub generate_obj {
    my $self = shift;
    return Genome::ProcessingProfile::SomaticVariation->create(@_);
}

sub get_required_params {
    my $class = shift;
    my $super = $class->SUPER::get_required_params;
    my @all = (@$super, @required_params);
    return \@all;
}

sub create_tiering_version {
    return 1;
}

1;

