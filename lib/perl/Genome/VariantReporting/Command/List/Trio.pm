package Genome::VariantReporting::Command::List::Trio;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Command::List::Trio {
    is => 'Genome::Object::Command::List',
    has => [
        subject_class_name => {
            is_constant => 1,
            value => 'Genome::VariantReporting::Process::Trio'
        },
        show => { default_value => 'id,status,created_by,metadata_directory' },
    ],
    doc => 'list trio processes',
};

1;
