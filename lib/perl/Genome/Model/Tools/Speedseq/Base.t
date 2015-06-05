#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 3;

my $pkg = 'Genome::Model::Tools::Speedseq::Base';
use_ok ($pkg);

my $version = 'test';
my $expected_path = '/gscmnt/sata849/info/speedseq_freeze/v1/speedseq/bin/speedseq';


my $path = $pkg->path_for_version($version);

ok($path eq $expected_path,'expected path for version');

$path = undef;
eval {
     $path = $pkg->path_for_version('foo');
};
ok(!$path,'no path for version');      