package Genome::Test::Factory::Model::RnaSeq;
use Genome::Test::Factory::Model;
@ISA = (Genome::Test::Factory::Model);

use strict;
use warnings;

use Genome;
use Genome::Test::Factory::Sample;
use Genome::Test::Factory::ProcessingProfile::RnaSeq;

our @required_params = qw(subject_id);

sub create_processing_profile_id {
    my $p = Genome::Test::Factory::ProcessingProfile::RnaSeq->setup_object();
    return $p->id;
}

sub generate_obj {
    my $self = shift;
    my $m = Genome::Model::RnaSeq->create(@_);
    return $m;
}

sub create_subject_id {
    my $subject = Genome::Test::Factory::Sample->setup_object();
    return $subject->id;
}

1;
