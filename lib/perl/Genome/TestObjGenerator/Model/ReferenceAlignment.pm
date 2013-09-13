package Genome::TestObjGenerator::Model::ReferenceAlignment;
use Genome::TestObjGenerator::Model;
@ISA = (Genome::TestObjGenerator::Model);

use strict;
use warnings;
use Genome;
use Genome::TestObjGenerator::ProcessingProfile::ReferenceAlignment;
use Genome::TestObjGenerator::Model::ReferenceSequence;
use Genome::TestObjGenerator::Build;

our @required_params = qw(reference_sequence_build subject_name subject_type);

sub generate_obj {
    my $self = shift;

    my $m = Genome::Model::ReferenceAlignment->create(@_);
    return $m;
}

sub create_processing_profile_id {
    my $p = Genome::TestObjGenerator::ProcessingProfile::ReferenceAlignment->setup_object();
    return $p->id;
}

sub create_reference_sequence_build {
    my $m = Genome::TestObjGenerator::Model::ReferenceSequence->setup_object;
    my $b = Genome::TestObjGenerator::Build->setup_object(model_id => $m->id);
    return $b;
}

sub create_subject_name {
    return Genome::TestObjGenerator::Util::generate_name("test_subject");
}

sub create_subject_type {
    return "sample_name";
}

1;
