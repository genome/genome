package Genome::DataSource::Oltp;

use strict;
use warnings;
use Genome;

class Genome::DataSource::Oltp {
    is => ['Genome::DataSource::OracleType', 'Genome::DataSource::CommonRDBMS'],
    has_classwide_constant => [
        server  => { default_value => 'gscprod' },
        login   => { default_value => 'gscuser' },
        auth    => { default_value => 'g_user' },
        owner   => { default_value => 'GSC' },
    ],
};

sub table_and_column_names_are_upper_case { 1; }

1;

