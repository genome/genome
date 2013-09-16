package Genome::Test::Factory::Model::ImportedAnnotation;
use Genome::Test::Factory::Model;
@ISA = (Genome::Test::Factory::Model);

use strict;
use warnings;
use Genome;
use Genome::Test::Factory::ProcessingProfile::ImportedAnnotation;
use Genome::Test::Factory::Taxon;

our @required_params = qw(subject);

sub generate_obj {
    my $self = shift;
    my $m = Genome::Model::ImportedAnnotation->create(@_);
    return $m;
}

sub create_processing_profile_id {
    my $p = Genome::Test::Factory::ProcessingProfile::ImportedAnnotation->setup_object();
    return $p->id;
}

sub create_subject {
    return Genome::Test::Factory::Taxon->setup_object();
}

1;
