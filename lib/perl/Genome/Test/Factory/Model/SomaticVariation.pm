package Genome::Test::Factory::Model::SomaticVariation;
use Genome::Test::Factory::Model;
@ISA = (Genome::Test::Factory::Model);

use strict;
use warnings;
use Genome;
use Genome::Test::Factory::ProcessingProfile::SomaticVariation;
use Genome::Test::Factory::Model::ImportedAnnotation;
use Genome::Test::Factory::Build;

our @required_params = qw(normal_model tumor_model annotation_build);

sub generate_obj {
    my $self = shift;
    my $m = Genome::Model::SomaticVariation->create(@_);
    return $m;
}

sub create_processing_profile_id {
    my $p = Genome::Test::Factory::ProcessingProfile::SomaticVariation->setup_object();
    return $p->id;
}

sub create_annotation_build {
    my $a = Genome::Test::Factory::Model::ImportedAnnotation->setup_object;
    my $annotation_build = Genome::Test::Factory::Build->setup_object(model_id => $a->id);
    return $annotation_build;
}

1;
