package Genome::Model::Build::ReferenceAlignment::View::ReferenceCoverage::Html;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::ReferenceAlignment::View::ReferenceCoverage::Html {
    is => 'UR::Object::View::Default::Html',
    has_constant => [
        perspective => { value => 'reference-coverage' },
    ],
};

sub _generate_content {
    my $self = shift;
    return 'Hello World';
}
