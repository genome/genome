package Genome::Model::Tools::RefCov::ROI::Region;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::RefCov::ROI::Region {
    is => ['Genome::Model::Tools::RefCov::ROI::RegionI'],
    has_optional => [
        name => { is => 'String', },
        chrom => { is => 'String', },
    ],
};


1;
