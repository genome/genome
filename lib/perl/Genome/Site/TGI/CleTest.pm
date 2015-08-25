package Genome::Site::TGI::CleTest;

use strict;
use warnings;
use Genome;

sub get_config {
    my $file = __FILE__.".yaml";
    my ($config, undef, undef) = Genome::Config::Parser->parse($file);
    return $config;
}

1;
