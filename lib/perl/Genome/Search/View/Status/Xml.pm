
package Genome::Search::View::Status::Xml;

use strict;
use warnings;

class Genome::Search::View::Status::Xml {
    is => 'UR::Object::View::Default::Xml',
    has_constant => [
        perspective => { value => 'status' },
        toolkit => { value => 'xml' },
        default_aspects => {
            value => ['interpreter', 'genome_path', 'ur_path']
        }
    ]
};

1;
