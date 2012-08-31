package Genome::View::Status::Html;

use strict;
use warnings;

use Genome;

class Genome::View::Status::Html {
    is => 'UR::Object::View::Default::Html',
    is_abstract => 1,
    has_constant => [
        perspective => { value => 'status' },
    ],
    doc => 'Concrete classes that want to be able to transform their Status::Xml views to html should inherit from this class.'
};

1;
