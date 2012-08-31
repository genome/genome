package Genome::Model::Build::Germline;

use strict;
use warnings;
use Genome;

class Genome::Model::Build::Germline {
    is => 'Genome::Model::Build',
    has => [        
        source_build => {
            is => 'Genome::Model::Build::ReferenceAlignment',
            id_by => 'source_build_id',
        },
        source_build_id => {
            is => 'Text',
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'source_id' ],
            is_many => 0,
            is_mutable => 1,
        },

    ],
};

1;

