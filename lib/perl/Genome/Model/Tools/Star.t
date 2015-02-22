#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;

my $pkg = 'Genome::Model::Tools::Star';
use_ok($pkg);

my $default_version = $pkg->default_star_version;
is($default_version, '2.3.1z1', 'default star version is correct');

my $default_path = $pkg->path_for_star_version($default_version);
ok($default_path, 'got default star tool path');

done_testing();
