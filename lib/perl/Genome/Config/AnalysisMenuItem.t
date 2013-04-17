#!/usr/bin/env genome-perl
use strict;
use warnings;

$ENV{UR_DBI_NO_COMMIT} = 1;

use Test::More;
use above "Genome";
use Carp::Always;

my $class = 'Genome::Config::AnalysisMenuItem';

use_ok($class);

my $menu_item = Genome::Config::AnalysisMenuItem->create(name => 'test');
ok($menu_item, 'it creates itself successfully');
isa_ok($menu_item, $class, 'it returns the correct object type');

done_testing();
