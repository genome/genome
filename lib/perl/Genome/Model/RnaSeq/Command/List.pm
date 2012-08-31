package Genome::Model::RnaSeq::Command::List;

use strict;
use warnings;

use Genome;
use Command;
use Data::Dumper;

class Genome::Model::RnaSeq::Command::List {
    is => 'UR::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1,
            value => 'Genome::Model::RnaSeq'
        },
        show => { default_value => 'id,name,processing_profile' },
    ],
    doc => 'list small-rna genome models',
};

1;

