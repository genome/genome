package Genome::DataSource::Dwrac;

use strict;
use warnings;
use Genome;

class Genome::DataSource::Dwrac {
    is => ['Genome::DataSource::OracleType'],
    has_classwide_constant => [
        server  => { default_value => 'dwrac' },
        login   => { default_value => 'gscuser' },
        auth    => { default_value => 'user_dw' },
        owner   => { default_value => 'GSC' },
    ],
};

sub table_and_column_names_are_upper_case { 1; }

1;

