package Genome::TestObjGenerator::ProcessingProfile::RnaSeq;
use Genome::TestObjGenerator::ProcessingProfile;
@ISA = (Genome::TestObjGenerator::ProcessingProfile);

use strict;
use warnings;
use Genome;

sub generate_obj {
    my $self = shift;
    my $p = Genome::ProcessingProfile::RnaSeq->create(@_);
    return $p;
}

1;

