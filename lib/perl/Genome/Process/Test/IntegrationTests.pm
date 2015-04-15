package Genome::Process::Test::IntegrationTests;

use strict;
use warnings;
use Genome;
use YAML;

sub get_process_test_configuration {
    my $process_name = shift;

    my ($configurations) = YAML::LoadFile(__FILE__ . '.YAML');
    return $configurations->{$process_name};
}

1;
