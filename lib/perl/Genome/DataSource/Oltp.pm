package Genome::DataSource::Oltp;

use strict;
use warnings;
use Genome;

class Genome::DataSource::Oltp {
    is => [ Genome::Config::get('ds_oltp_type'), 'Genome::DataSource::CommonRDBMS' ],
    has_classwide_constant => [
        server  => { default_value => Genome::Config::get('ds_oltp_server') },
        login   => { default_value => Genome::Config::get('ds_oltp_login') },
        auth    => { default_value => Genome::Config::get('ds_oltp_auth') },
        owner   => { default_value => Genome::Config::get('ds_oltp_owner') },
    ],
};

sub table_and_column_names_are_upper_case { 0; }

1;

