#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::GenotypeMicroarray::Filter::ByGcScore') or die;

my $filter = Genome::Model::GenotypeMicroarray::Filter::ByGcScore->create();
ok(!$filter, 'create failed w/o min');
$filter = Genome::Model::GenotypeMicroarray::Filter::ByGcScore->create(min => -0.1);
ok(!$filter, 'create failed w/ invalid min');

$filter = Genome::Model::GenotypeMicroarray::Filter::ByGcScore->create(min => 0.7);
ok($filter, 'create filter');
ok($filter->filter({gc_score => undef}), 'did not filter undef');
ok($filter->filter({gc_score => -1}), 'did not filter -1');
ok($filter->filter({gc_score => .71}), 'did not filter .71');
ok($filter->filter({gc_score => .70}), 'did not filter .70');
ok(!$filter->filter({gc_score => .69}), 'filtered .69');

done_testing();
exit;

