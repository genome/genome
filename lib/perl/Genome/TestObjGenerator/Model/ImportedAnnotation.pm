package Genome::TestObjGenerator::Model::ImportedAnnotation;
use Genome::TestObjGenerator::Model;
@ISA = (Genome::TestObjGenerator::Model);

use strict;
use warnings;
use Genome;
use Genome::TestObjGenerator::ProcessingProfile::ImportedAnnotation;
use Genome::TestObjGenerator::Taxon;

my @required_params = ("subject");

sub generate_obj {
    my $self = shift;
    my $m = Genome::Model::ImportedAnnotation->create(@_);
    return $m;
}

sub create_processing_profile_id {
    my $p = Genome::TestObjGenerator::ProcessingProfile::ImportedAnnotation->setup_object();
    return $p->id;
}

sub create_subject {
    return Genome::TestObjGenerator::Taxon->setup_object();
}

sub get_required_params {
    my $class = shift;
    my $super = $class->SUPER::get_required_params;
    my @all = (@$super, @required_params);
    return \@all;
}

1;

