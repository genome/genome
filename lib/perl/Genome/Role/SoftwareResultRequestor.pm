package Genome::Role::SoftwareResultRequestor;

use strict;
use warnings;

use Genome;
use UR::Role;

# This role is intended to be mixed into classes that can act as the requestor
# of a software result. A "requestor" is the entity that required the creation
# of the result
role Genome::Role::SoftwareResultRequestor { };

1;
