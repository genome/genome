package Genome::DataSource::Oltp;

use strict;
use warnings;
use Genome;

class Genome::DataSource::Oltp {
    is => [$ENV{GENOME_DS_OLTP_TYPE}, 'Genome::DataSource::CommonRDBMS'],
    has_classwide_constant => [
        server  => { default_value => $ENV{GENOME_DS_OLTP_SERVER} },
        login   => { default_value => $ENV{GENOME_DS_OLTP_LOGIN} },
        auth    => { default_value => $ENV{GENOME_DS_OLTP_AUTH}  },
        owner   => { default_value => $ENV{GENOME_DS_OLTP_OWNER} },
    ],
};

sub table_and_column_names_are_upper_case { 1; }

1;

