package Genome::Test::Factory::Model::ClinSeq;
use Genome::Test::Factory::Model;
@ISA = (Genome::Test::Factory::Model);

use strict;
use warnings;
use Genome;
use Genome::Test::Factory::ProcessingProfile::ClinSeq;

sub generate_obj {
    my $self = shift;
    my $m = Genome::Model::ClinSeq->create(@_);
    return $m;
}

sub create_processing_profile_id {
    my $p = Genome::Test::Factory::ProcessingProfile::ClinSeq->setup_object();
    return $p->id;
}

1;
