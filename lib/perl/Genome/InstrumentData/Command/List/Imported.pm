package Genome::InstrumentData::Command::List::Imported;

use strict;
use warnings;

use Genome;
use Command;
use Data::Dumper;

class Genome::InstrumentData::Command::List::Imported {
    is => 'Genome::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1, 
            value => 'Genome::InstrumentData::Imported' 
        },
        show => { default_value => 'id,sample_name,sequencing_platform,import_format' },
    ],
    doc => 'list imported instrument data available for analysis',
};

1;
