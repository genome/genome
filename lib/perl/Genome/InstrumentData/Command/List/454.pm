package Genome::InstrumentData::Command::List::454;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::List::454 {
    is => 'Genome::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1,
            value => 'Genome::InstrumentData::454'
        },
        show => { default_value => 'id,run_name,region_number,index_sequence,sample_name' }
    ],
    doc => 'list 454 regions available for analysis',
};

sub _base_filter {
    'sequencing_platform=454'
}

1;

