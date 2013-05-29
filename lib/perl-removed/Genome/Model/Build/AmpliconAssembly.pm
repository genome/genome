package Genome::Model::Build::AmpliconAssembly;

use strict;
use warnings;

use Genome;

use Carp 'confess';

class Genome::Model::Build::AmpliconAssembly {
    is => 'Genome::Model::Build',
};

sub create {
    die __PACKAGE__ . ' is deprecated.';
}

1;
