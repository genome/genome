package Genome::Model::ReferenceSequence::Command::Converter::List;

use strict;
use warnings;

use Genome;

class Genome::Model::ReferenceSequence::Command::Converter::List {
    is => 'UR::Object::Command::List',
    has => [
        subject_class_name => {
            is_constant => 1,
            value => 'Genome::Model::Build::ReferenceSequence::Converter'
        },
        show => { default_value => 'id,source_reference_build,destination_reference_build,algorithm,resource'},
    ],
    doc => 'list reference-sequence converters',
};

1;
