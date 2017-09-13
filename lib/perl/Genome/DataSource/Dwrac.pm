package Genome::DataSource::Dwrac;

use strict;
use warnings;
use Genome;

class Genome::DataSource::Dwrac {
    is => [ Genome::Config::get('ds_dwrac_type'), 'Genome::DataSource::CommonRDBMS' ],
    has_classwide_constant => [
        server  => { default_value => Genome::Config::get('ds_dwrac_server') },
        login   => { default_value => Genome::Config::get('ds_dwrac_login') },
        auth    => { default_value => Genome::Config::get('ds_dwrac_auth') },
        owner   => { default_value => Genome::Config::get('ds_dwrac_owner') },
    ],
};

sub table_and_column_names_are_upper_case { 0; }

1;

