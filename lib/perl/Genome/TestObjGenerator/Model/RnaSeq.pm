package Genome::TestObjGenerator::Model::RnaSeq;
use Genome::TestObjGenerator::Model;
@ISA = (Genome::TestObjGenerator::Model);

use strict;
use warnings;
use Genome;
use Genome::TestObjGenerator::ProcessingProfile::RnaSeq;

sub create_processing_profile_id {
    my $p = Genome::TestObjGenerator::ProcessingProfile::RnaSeq->setup_object();
    return $p->id;
}

sub generate_obj {
    my $self = shift;
    my $m = Genome::Model::RnaSeq->create(@_);
    return $m;
}

1;

