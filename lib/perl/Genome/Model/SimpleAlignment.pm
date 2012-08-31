package Genome::Model::SimpleAlignment;

use strict;
use warnings;

use Genome;

class Genome::Model::SimpleAlignment {
    is  => 'Genome::ModelDeprecated',
    has => [
       reference_sequence_name => { via => 'processing_profile'},
    ],
    
};


1
;
