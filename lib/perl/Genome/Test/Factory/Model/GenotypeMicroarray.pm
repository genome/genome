package Genome::Test::Factory::Model::GenotypeMicroarray;
use Genome::Test::Factory::Model;
@ISA = (Genome::Test::Factory::Model);

use strict;
use warnings;

use Genome;
use Genome::Test::Factory::ProcessingProfile::GenotypeMicroarray;
use Genome::Test::Factory::Model::ImportedVariationList;
use Genome::Test::Factory::Build;
use Genome::Test::Factory::Sample;

our @required_params = qw(dbsnp_build subject_id);

sub generate_obj {
    my $self = shift;

    my $m = Genome::Model::GenotypeMicroarray->create(@_);
    return $m;
}

sub create_processing_profile_id {
    my $p = Genome::Test::Factory::ProcessingProfile::GenotypeMicroarray->setup_object();
    return $p->id;
}

sub create_dbsnp_build {
    my $m = Genome::Test::Factory::Model::ImportedVariationList->setup_object;
    my $b = Genome::Test::Factory::Build->setup_object(model_id => $m->id);
    return $b;
}

sub create_subject_id {
    my $subject = Genome::Test::Factory::Sample->setup_object();
    return $subject->id;
}

1;
