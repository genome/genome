use strict;
use warnings;

package Genome::DataSource::OldGMSchema;
use Genome;

class Genome::DataSource::OldGMSchema {
    is => ['UR::DataSource::Oracle'],
    has_constant => [
        server => { default_value => 'dwrac' },
        login => { default_value => 'mguser' },
        auth => { default_value => 'mguser_prd' },
        owner => { default_value => 'MG' },
    ],
};

sub table_and_column_names_are_upper_case { 1; }

1;
