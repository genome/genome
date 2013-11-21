package Genome::DataSource::Dwrac;

use strict;
use warnings;
use Genome;

class Genome::DataSource::Dwrac {
    is => [$ENV{GENOME_DS_DWRAC_TYPE}, 'Genome::DataSource::CommonRDBMS'],
    has_classwide_constant => [
        server  => { default_value => $ENV{GENOME_DS_DWRAC_SERVER} },
        login   => { default_value => $ENV{GENOME_DS_DWRAC_LOGIN} },
        auth    => { default_value => $ENV{GENOME_DS_DWRAC_AUTH} },
        owner   => { default_value => $ENV{GENOME_DS_DWRAC_OWNER} },
    ],
};

sub table_and_column_names_are_upper_case { 1; }

1;

