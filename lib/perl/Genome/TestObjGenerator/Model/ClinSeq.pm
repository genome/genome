package Genome::TestObjGenerator::Model::ClinSeq;
use Genome::TestObjGenerator::Model;
@ISA = (Genome::TestObjGenerator::Model);

use strict;
use warnings;
use Genome;
use Genome::TestObjGenerator::ProcessingProfile::ClinSeq;

sub generate_obj {
    my $self = shift;
    my $m = Genome::Model::ClinSeq->create(@_);
    return $m;
}

sub create_processing_profile_id {
    my $p = Genome::TestObjGenerator::ProcessingProfile::ClinSeq->setup_object();
    return $p->id;
}

1;
