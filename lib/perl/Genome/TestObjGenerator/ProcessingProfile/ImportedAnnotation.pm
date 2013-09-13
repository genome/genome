package Genome::TestObjGenerator::ProcessingProfile::ImportedAnnotation;
use Genome::TestObjGenerator::ProcessingProfile;
@ISA = (Genome::TestObjGenerator::ProcessingProfile);

use strict;
use warnings;

my @required_params = ("annotation_source");

sub get_required_params {
    my $class = shift;
    my $super = $class->SUPER::get_required_params;
    my @all = (@$super, @required_params);
    return \@all;
}

sub create_annotation_source {
    return Genome::TestObjGenerator::Util::generate_name("test_annotation_source");
}

1;
