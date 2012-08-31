#!/usr/bin/env genome-perl

use strict;
use warnings;
use above "Genome"; 
use Test::More tests => 4;

my $event_id = 90005615;
my $event = Genome::Model::Event->get($event_id);
ok($event, "got an event");

my $m = $event->add_metric(
    name => 'myname',
    value => 'myvalue',
);
ok($m, "made a metric");
is($m->name, 'myname', 'name is correct');
is($m->value, 'myvalue', 'value is correct');


