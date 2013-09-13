package Genome::TestObjGenerator::Model::SomaticVariation;
use Genome::TestObjGenerator::Model;
@ISA = (Genome::TestObjGenerator::Model);

use strict;
use warnings;
use Genome;
use Genome::TestObjGenerator::ProcessingProfile::SomaticVariation;
use Genome::TestObjGenerator::Model::ImportedAnnotation;
use Genome::TestObjGenerator::Build;

our @required_params = qw(normal_model tumor_model annotation_build);

sub generate_obj {
    my $self = shift;
    my $m = Genome::Model::SomaticVariation->create(@_);
    return $m;
}

sub create_processing_profile_id {
    my $p = Genome::TestObjGenerator::ProcessingProfile::SomaticVariation->setup_object();
    return $p->id;
}

sub create_annotation_build {
    my $a = Genome::TestObjGenerator::Model::ImportedAnnotation->setup_object;
    my $annotation_build = Genome::TestObjGenerator::Build->setup_object(model_id => $a->id);
    return $annotation_build;
}

1;
