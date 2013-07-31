package Genome::TestObjGenerator::ProcessingProfile::ImportedVariationList;
use Genome::TestObjGenerator::ProcessingProfile;
@ISA = (Genome::TestObjGenerator::ProcessingProfile);

use strict;
use warnings;
use Genome;

sub generate_obj {
    my $self = shift;

    my $p = Genome::ProcessingProfile::ImportedVariationList->create(@_);
    return $p;
}

1;

