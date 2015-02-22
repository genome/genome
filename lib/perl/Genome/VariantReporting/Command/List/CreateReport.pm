package Genome::VariantReporting::Command::List::CreateReport;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Command::List::CreateReport {
    is => 'Genome::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1,
            value => 'Genome::VariantReporting::Process::CreateReport'
        },
        show => { default_value => 'id,status,created_by,metadata_directory' },
    ],
    doc => 'list create-report processes',
};

1;
