package Genome::Role::SoftwareResultSponsor;

use strict;
use warnings;

use Genome;
use UR::Role;

# This role is intended to be mixed into classes that can act as the sponsor
# of a software result. A "sponsor" is the entity that is accountable for
# the resource usage of the result.
role Genome::Role::SoftwareResultSponsor { };

1;
