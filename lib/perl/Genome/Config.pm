package Genome::Config;

use strict;
use warnings;

sub get {
    my $config_key = shift;
    my $env_key = 'GENOME_' . uc($config_key);
    return $ENV{$env_key};
}

1;
