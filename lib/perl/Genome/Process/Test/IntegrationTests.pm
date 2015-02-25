package Genome::Process::Test::IntegrationTests;

use strict;
use warnings;
use Genome;
use YAML;

sub process_yaml {
    my $process_name = shift;

    my ($yaml) = YAML::LoadFile(__FILE__ . '.YAML');
    return $yaml->{$process_name};
}

1;
