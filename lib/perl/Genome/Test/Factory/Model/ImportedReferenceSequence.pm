package Genome::Test::Factory::Model::ImportedReferenceSequence;
use Genome::Test::Factory::Model;
@ISA = (Genome::Test::Factory::Model);

use strict;
use warnings;
use Genome;
use Genome::Test::Factory::ProcessingProfile::ImportedReferenceSequence;

our @required_params = qw(subject);

sub generate_obj {
    my $self = shift;

    my $m = Genome::Model::ImportedReferenceSequence->create(@_);
    return $m;
}

sub create_processing_profile_id {
    my $p = Genome::Test::Factory::ProcessingProfile::ImportedReferenceSequence->setup_object();
    return $p->id;
}

sub create_subject {
    return Genome::Sample->create(name => Genome::Test::Factory::Util::generate_name("test_subject"));
}

1;
