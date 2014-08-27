package Genome::SoftwareResult::Input;

use strict;
use warnings;

use Genome;

class Genome::SoftwareResult::Input {
    is => 'Genome::SoftwareResult::PIBase',
    table_name => 'result.input',
    type_name => 'software result input',
    id_by => [
        software_result => {
            is => 'Genome::SoftwareResult',
            id_by => 'software_result_id',
        },
        name => {
            is => 'Text',
            len => 255,
            column_name => 'input_name',
        },
    ],
    has => [
        value_id => {
            is => 'Text',
            len => 1000,
            column_name => 'input_value',
        },
        value_obj => {
            is => 'UR::Object',
            id_by => 'value_id',
            id_class_by => 'value_class_name',
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

sub input_name {
    Carp::confess("using input_name!");
}

sub input_value {
    Carp::confess("using input_value!");
}

1;
