package Genome::SoftwareResult::Requestor;

use strict;
use warnings;

use Genome;

class Genome::SoftwareResult::Requestor {
    is => 'UR::Object',
    is_abstract => 1,
    doc => 'This "role" is intended to be mixed into classes that can act as the requestor of a software result. A "requestor" is the entity that required the creation of the result',
};


1;
