package Genome::View::Resource::Html;

use strict;
use warnings;

use Genome;

class Genome::View::Resource::Html {
    is => 'UR::Object::View::Default::Html',
    is_abstract => 1,
    has_constant => [
        perspective => 'resource',
    ],
    doc => 'Placeholder class so that the rest app can resolve this view to find static resources relative to it'
};

1;
