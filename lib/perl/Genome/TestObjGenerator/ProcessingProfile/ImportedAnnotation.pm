package Genome::TestObjGenerator::ProcessingProfile::ImportedAnnotation;
use Genome::TestObjGenerator::ProcessingProfile;
@ISA = (Genome::TestObjGenerator::ProcessingProfile);

use strict;
use warnings;

our @required_params = qw(annotation_source);

sub create_annotation_source {
    return Genome::TestObjGenerator::Util::generate_name("test_annotation_source");
}

1;
