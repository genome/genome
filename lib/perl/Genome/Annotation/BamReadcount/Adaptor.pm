package Genome::Annotation::BamReadcount::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::BamReadcount::Adaptor {
    is => 'Genome::Annotation::AdaptorBase',

    has_planned_output => [
        version => { is  => 'Version', },
        minimum_mapping_quality => { is => 'Integer', },
        minimum_base_quality => { is => 'Integer', },
        max_count => { is  => 'Integer', },
        per_library => { is  => 'Bool', },
        insertion_centric => { is  => 'Bool', },
    ],
};


1;
