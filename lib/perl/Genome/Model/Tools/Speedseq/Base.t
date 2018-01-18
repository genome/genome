#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 5;

my $pkg = 'Genome::Model::Tools::Speedseq::Base';
use_ok ($pkg);

my $version = '0.1.0-gms';
my $expected_path = '/gscmnt/ams1102/info/speedseq_freeze/v3/speedseq/bin/speedseq';


my $path = $pkg->path_for_version($version);

ok($path eq $expected_path,'expected path for version');

$path = undef;
eval {
     $path = $pkg->path_for_version('foo');
};
ok(!$path,'no path for version');

is(Genome::Model::Tools::Speedseq::Base::_format_tool_arg('Boolean', 'test', 0), '', "False boolean argument formats as an empty string");
is(Genome::Model::Tools::Speedseq::Base::_format_tool_arg('Boolean', 'test', 1), '-test', "True boolean argument formats as a dash with the arg name");
