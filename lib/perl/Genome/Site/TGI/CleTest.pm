package Genome::Site::TGI::CleTest;

use strict;
use warnings;
use Genome;
use YAML;

sub get_config {
    my $file = __FILE__.".yaml";
    my ($config, undef, undef) = YAML::LoadFile($file);
    return $config;
}

1;
