#!/usr/bin/env genome-perl

use Test::More tests => 3;
use Data::Dumper;

use GSCApp; #use App; # This should work, but App won't stand-alone. FIXME
App->init; 

$o1 = Bio::Gene->get("foo_gene"); 
ok($o1, "got a nonsense gene");

$o2 = Bio::Gene->get("foo_gene");
ok($o2, "got the same nonsense gene on the second try");

ok($o1 eq $o2, "both are the same reference");
