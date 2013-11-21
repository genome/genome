#$Id$

package Genome::Model::GenePrediction::Command::Pap::ListChunker;

use strict;
use warnings;

use Workflow;

use List::MoreUtils qw/part/;
use English;
use Carp;


class Genome::Model::GenePrediction::Command::Pap::ListChunker {
    is  => ['Command::V1'],
    has => [
        list  => { is => 'ARRAY', doc => 'list to split up'                             },
        chunk_number  => { is => 'SCALAR', doc => 'number of sequences per output file'         },
        list_of_lists  => { is => 'ARRAY', doc => 'list of lists/array of arrays that are split out',
                            is_optional => 1,         },
    ],
};

1;
