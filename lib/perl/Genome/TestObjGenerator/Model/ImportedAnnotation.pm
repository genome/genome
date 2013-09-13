package Genome::TestObjGenerator::Model::ImportedAnnotation;
use Genome::TestObjGenerator::Model;
@ISA = (Genome::TestObjGenerator::Model);

use strict;
use warnings;
use Genome;
use Genome::TestObjGenerator::ProcessingProfile::ImportedAnnotation;
use Genome::TestObjGenerator::Taxon;

our @required_params = qw(subject);

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

1;
