#!/usr/bin/env genome-perl

use strict;
use warnings;
use above 'Genome';
use Test::More tests => 1;

class Genome::Model::TestPipeLine {
    is => 'Genome::Model'
};

my $dir = $INC{"Genome.pm"};
use File::Basename;
$dir = dirname($dir);
chdir $dir;

my $exit_code = system "genome processing-profile list small-rna 2>/dev/null 1>/dev/null";
$exit_code /= 256;
is($exit_code, 0, "list for the small-rna processing profiles works");


