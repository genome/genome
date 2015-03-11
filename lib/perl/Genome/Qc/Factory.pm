package Genome::Qc::Factory;

use strict;
use warnings;
use Genome;

sub get_config {
    return Genome::Qc::Config->create(@_);
}

1;

