package Genome::InstrumentData::AlignmentResult::Command::List;
use strict;
use warnings;
use Genome;

class Genome::InstrumentData::AlignmentResult::Command::List {
    is => 'UR::Object::Command::List',
    has => [
        subject_class_name => { is_constant => 1, value => 'Genome::InstrumentData::AlignmentResult' },
        filter => { 
            default_value => '' 
        },
        show => {
            default_value => 'id,instrument_data_id,reference_name,aligner,trimmer,filter,module_version,test_name,output_dir'
        },
    ],
    doc => 'list alignment data sets'
};

1;

