#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More tests => 1;
use File::Path;


BEGIN
{
    use_ok ('Genome::Model::Tools::Fasta::Synch');
}

1;
