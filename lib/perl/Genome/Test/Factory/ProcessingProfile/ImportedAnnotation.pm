package Genome::Test::Factory::ProcessingProfile::ImportedAnnotation;
use Genome::Test::Factory::ProcessingProfile;
@ISA = (Genome::Test::Factory::ProcessingProfile);

use strict;
use warnings;

our @required_params = qw(annotation_source);

sub create_annotation_source {
    return Genome::Test::Factory::Util::generate_name("test_annotation_source");
}

1;
