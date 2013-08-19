package Genome::TestObjGenerator::ProcessingProfile::SomaticValidation;
use Genome::TestObjGenerator::ProcessingProfile;
@ISA = (Genome::TestObjGenerator::ProcessingProfile);

use strict;
use warnings;
use Genome;

sub generate_obj {
    my $self = shift;
    return Genome::ProcessingProfile::SomaticValidation->create(@_);
}

1;

