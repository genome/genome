package Genome::TestObjGenerator::Model::SomaticVariation;
use Genome::TestObjGenerator::Model;
@ISA = (Genome::TestObjGenerator::Model);

use strict;
use warnings;
use Genome;
use Genome::TestObjGenerator::ProcessingProfile::SomaticVariation;
use Genome::TestObjGenerator::Model::ImportedAnnotation;
use Genome::TestObjGenerator::Build;

my @required_params = ("normal_model", "tumor_model", "annotation_build", "subject_name");

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

sub get_required_params {
    my $class = shift;
    my $super = $class->SUPER::get_required_params;
    my @all = (@$super, @required_params);
    return \@all;
}

1;

