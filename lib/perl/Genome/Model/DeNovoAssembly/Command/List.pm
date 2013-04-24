package Genome::Model::DeNovoAssembly::Command::List;

use strict;
use warnings;

use Genome;
use Command; 
use Data::Dumper;

class Genome::Model::DeNovoAssembly::Command::List {
    is => 'Genome::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1, 
            value => 'Genome::Model::DeNovoAssembly' 
        },
        show => { default_value => 'id,name,subject,read_processor,assembler_name,assembler_version,assembler_params' },
    ],
    doc => 'list assembly genome models',
};

1;

