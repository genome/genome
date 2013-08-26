package Genome::DataSource::Dgidb;

use strict;
use warnings;
use Genome;

class Genome::DataSource::Dgidb {
    is => 'UR::DataSource::Pg',
    has_constant => [
        server => { default_value => 'dbname=genome;host=postgres' },
        login => { default_value => 'genome' },
        auth => { default_value => 'TGI_pg_1' },
        owner => { default_value => 'public' },
    ],
};

sub table_and_column_names_are_upper_case { 0 }

1;
