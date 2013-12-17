package Genome::DataSource::Dgidb;

use strict;
use warnings;
use Genome;

class Genome::DataSource::Dgidb {
    is => $ENV{GENOME_DS_DGIDB_TYPE},
    has_constant => [
        server  => { default_value => $ENV{GENOME_DS_DGIDB_SERVER} },
        login   => { default_value => $ENV{GENOME_DS_DGIDB_LOGIN} },
        auth    => { default_value => $ENV{GENOME_DS_DGIDB_AUTH} },
        owner   => { default_value => $ENV{GENOME_DS_DGIDB_OWNER} },
    ],
};

sub table_and_column_names_are_upper_case { 0 }

1;
