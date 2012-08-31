package Genome::DataSource::PGTest;

use strict;
use warnings;
use Genome;
use Carp;

class Genome::DataSource::PGTest {
    is => 'Genome::DataSource::Main',
    has_constant => [
        server => { default_value => 'dbname=genome;host=gms-postgres' },
        login => { default_value => 'genome' },
        auth => { default_value => 'TGIlab' },
        owner => { default_value => 'public' },
    ],
};

sub _ds_tag {
    'Genome::DataSource::PGTest';
}

1;
