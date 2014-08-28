package Genome::SoftwareResult::Param;

use strict;
use warnings;

use Genome;

class Genome::SoftwareResult::Param {
    is => 'Genome::SoftwareResult::PIBase',
    table_name => 'result.param',
    has => [
        value_id => {
            is => 'Text',
            len => 1000,
            column_name => 'param_value',
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

sub param_name {
    Carp::confess("using param_name!");
}

sub param_value {
    Carp::confess("using param_value!");
}

1;
