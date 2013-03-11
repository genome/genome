#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use File::Remove qw/ remove /;

use Test::More tests => 1;


BEGIN {
    use_ok('Genome::Model::Tools::Ber::AmgapBerProtName');
}

unless( -d "/tmp/disk/BER_TEST")
{
    mkdir("/tmp/disk/BER_TEST");
}

unless( -l "/tmp/disk/analysis")
{
    symlink("/gscmnt/temp110/info/annotation/wnash/BER_TEST",
            "/tmp/disk/BER_TEST");
}
