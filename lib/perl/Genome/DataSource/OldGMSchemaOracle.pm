package Genome::DataSource::OldGMSchemaOracle;

use strict;
use warnings;
use Genome;

class Genome::DataSource::OldGMSchemaOracle {
    is => ['Genome::DataSource::OracleType', 'Genome::DataSource::CommonRDBMS'],
    has_classwide_constant => [
        server  => { default_value => 'dwrac' },
        login   => { default_value => 'mguser' },
        auth    => { default_value => 'mguser_prd' },
        owner   => { default_value => 'MG' },
    ],
};

sub table_and_column_names_are_upper_case { 1; }

1;

