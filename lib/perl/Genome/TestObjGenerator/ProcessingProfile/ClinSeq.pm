package Genome::TestObjGenerator::ProcessingProfile::ClinSeq;
use Genome::TestObjGenerator::ProcessingProfile;
@ISA = (Genome::TestObjGenerator::ProcessingProfile);

use strict;
use warnings;
use Genome;

sub generate_obj {
    my $self = shift;

    my $p = Genome::ProcessingProfile::ClinSeq->create(@_);
    return $p
}

1;

