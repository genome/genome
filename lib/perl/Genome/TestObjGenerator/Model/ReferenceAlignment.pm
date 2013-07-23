package Genome::TestObjGenerator::Model::ReferenceAlignment;
use Genome::TestObjGenerator::Model;
@ISA = (Genome::TestObjGenerator::Model);

use strict;
use warnings;
use Genome;
use Genome::TestObjGenerator::ProcessingProfile::ReferenceAlignment;

my @required_params = ("reference_sequence_build");

sub generate_obj {
    my $self = shift;

    my $m = Genome::Model::ReferenceAlignment->create(@_);
    return $m;
}

sub create_processing_profile_id {
    my $p = Genome::TestObjGenerator::ProcessingProfile::ReferenceAlignment->setup_object();
    return $p->id;
}

sub create_reference_sequence_build {
#TODO Must fix before pushing!
    my $b = Genome::Model::Build::ReferenceSequence->get_by_name('NCBI-human-build36');
    return $b;
}

sub get_required_params {
    my $class = shift;
    my $super = $class->SUPER::get_required_params;
    my @all = (@$super, @required_params);
    return \@all;
}

1;

