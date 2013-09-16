package Genome::Test::Factory::Model::ReferenceAlignment;
use Genome::Test::Factory::Model;
@ISA = (Genome::Test::Factory::Model);

use strict;
use warnings;
use Genome;
use Genome::Test::Factory::ProcessingProfile::ReferenceAlignment;
use Genome::Test::Factory::Model::ReferenceSequence;
use Genome::Test::Factory::Build;

our @required_params = qw(reference_sequence_build subject_name subject_type);

sub generate_obj {
    my $self = shift;

    my $m = Genome::Model::ReferenceAlignment->create(@_);
    return $m;
}

sub create_processing_profile_id {
    my $p = Genome::Test::Factory::ProcessingProfile::ReferenceAlignment->setup_object();
    return $p->id;
}

sub create_reference_sequence_build {
    my $m = Genome::Test::Factory::Model::ReferenceSequence->setup_object;
    my $b = Genome::Test::Factory::Build->setup_object(model_id => $m->id);
    return $b;
}

sub create_subject_name {
    return Genome::Test::Factory::Util::generate_name("test_subject");
}

sub create_subject_type {
    return "sample_name";
}

1;
