package Genome::TestObjGenerator::ProcessingProfile::DifferentialExpression;
use Genome::TestObjGenerator::ProcessingProfile;
@ISA = (Genome::TestObjGenerator::ProcessingProfile);

use strict;
use warnings;
use Genome;

sub generate_obj {
    my $self = shift;
    my $p = Genome::ProcessingProfile::DifferentialExpression->create(@_);
    return $p;
}

1;

