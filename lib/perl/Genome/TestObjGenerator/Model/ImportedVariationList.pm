package Genome::TestObjGenerator::Model::ImportedVariationList;
use Genome::TestObjGenerator::Model;
@ISA = (Genome::TestObjGenerator::Model);

use strict;
use warnings;
use Genome;

use Genome::TestObjGenerator::ProcessingProfile::ImportedVariationList;

our @required_params = qw(subject);

sub generate_obj {
    my $self = shift;

    my $m = Genome::Model::ImportedVariationList->create(@_);
    return $m;
}

sub create_processing_profile_id {
    my $p = Genome::TestObjGenerator::ProcessingProfile::ImportedVariationList->setup_object();
    return $p->id;
}

sub create_subject {
    return Genome::Sample->create(name => Genome::TestObjGenerator::Util::generate_name("test_subject"));
}

1;
