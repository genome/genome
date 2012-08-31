package Genome::Model::Tools::Gtf::Base;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Gtf::Base {
    is  => 'Command::V2',
    is_abstract => 1,
    has => [
        input_gtf_file => {
            is => 'Text',
            doc => 'The input gtf format file to operate on.',
        },
    ],
};


1;
