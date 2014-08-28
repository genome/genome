package Genome::SoftwareResult::Input;

use strict;
use warnings;

use Genome;

class Genome::SoftwareResult::Input {
    is => 'Genome::SoftwareResult::PIBase',
    table_name => 'result.input',
    has => [
        value_id => {
            is => 'Text',
            len => 1000,
            column_name => 'input_value',
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
