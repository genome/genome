package Genome::DataSource::Dgidb;

use strict;
use warnings;
use Genome;

class Genome::DataSource::Dgidb {
    is => Genome::Config::get('ds_dgidb_type'),
    has_constant => [
        server  => { default_value => Genome::Config::get('ds_dgidb_server') },
        login   => { default_value => Genome::Config::get('ds_dgidb_login') },
        auth    => { default_value => Genome::Config::get('ds_dgidb_auth') },
        owner   => { default_value => Genome::Config::get('ds_dgidb_owner') },
    ],
};

sub table_and_column_names_are_upper_case { 0 }

1;
