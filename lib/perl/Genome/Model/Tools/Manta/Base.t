#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 3;
use Genome::Utility::Test qw(compare_ok);

my $pkg = 'Genome::Model::Tools::Manta::Base';
use_ok ($pkg);

my $version = '0.29.6';
my $expected_path = '/gscmnt/gc13001/info/model_data/jwalker_scratch/src/manta-0.29.6.centos5_x86_64/bin';

my $path = $pkg->path_for_version($version);

is($path,$expected_path,'expected path for version');

$path = undef;
eval {
     $path = $pkg->path_for_version('foo');
};
ok(!$path,'no path for version'); 