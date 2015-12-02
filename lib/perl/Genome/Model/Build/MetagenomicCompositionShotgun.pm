package Genome::Model::Build::MetagenomicCompositionShotgun;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::MetagenomicCompositionShotgun{
    is => 'Genome::Model::Build',
};
# This build type has 'from_build_links' 

sub create {
    die __PACKAGE__ . ' is deprecated.';
}

1;

