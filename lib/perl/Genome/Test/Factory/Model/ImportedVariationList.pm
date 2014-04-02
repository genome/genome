package Genome::Test::Factory::Model::ImportedVariationList;
use Genome::Test::Factory::Model;
@ISA = (Genome::Test::Factory::Model);

use strict;
use warnings;
use Genome;

use Genome::Test::Factory::ProcessingProfile::ImportedVariationList;
use Genome::Test::Factory::Model::ReferenceSequence;
use Genome::Test::Factory::Build;

our @required_params = qw(subject reference);

sub generate_obj {
    my $self = shift;

    my $m = Genome::Model::ImportedVariationList->create(@_);
    return $m;
}

sub create_processing_profile_id {
    my $p = Genome::Test::Factory::ProcessingProfile::ImportedVariationList->setup_object();
    return $p->id;
}

sub create_subject {
    return Genome::Sample->create(name => Genome::Test::Factory::Util::generate_name("test_subject"));
}

sub create_reference {
    my $m = Genome::Test::Factory::Model::ReferenceSequence->setup_object();
    my $b = Genome::Test::Factory::Build->setup_object(model_id => $m->id);
    return $b;
}

1;
