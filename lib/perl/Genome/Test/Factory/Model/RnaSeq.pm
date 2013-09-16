package Genome::Test::Factory::Model::RnaSeq;
use Genome::Test::Factory::Model;
@ISA = (Genome::Test::Factory::Model);

use strict;
use warnings;
use Genome;
use Genome::Test::Factory::ProcessingProfile::RnaSeq;

our @required_params = qw(subject_name subject_type);

sub create_processing_profile_id {
    my $p = Genome::Test::Factory::ProcessingProfile::RnaSeq->setup_object();
    return $p->id;
}

sub generate_obj {
    my $self = shift;
    my $m = Genome::Model::RnaSeq->create(@_);
    return $m;
}

sub create_subject_name {
    return Genome::Test::Factory::Util::generate_name("test_subject");
}

sub create_subject_type {
    return "sample_name";
}

1;
