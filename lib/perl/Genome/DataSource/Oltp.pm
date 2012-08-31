package Genome::DataSource::Oltp;

use strict;
use warnings;
use Genome;

class Genome::DataSource::Oltp {
    is => 'UR::DataSource::Oracle',
    type_name => 'genome datasource oltp',
};

sub server {
    'gscprod';
}

sub login {
    'gscuser';
}

sub auth {
    'g_user';
}

sub owner {
    'gsc';
}


1;

