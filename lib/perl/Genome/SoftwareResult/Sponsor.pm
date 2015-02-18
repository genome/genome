package Genome::SoftwareResult::Sponsor;

use strict;
use warnings;

use Genome;

class Genome::SoftwareResult::Sponsor {
    is => 'UR::Object',
    is_abstract => 1,
    doc => 'This "role" is intended to be mixed into classes that can act as the sponsor of a software result. A "sponsor" is the entity that is accountable for the resource usage of the result.',
};


1;
