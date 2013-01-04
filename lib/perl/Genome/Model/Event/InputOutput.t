#!/usr/bin/env genome-perl

use strict;
use warnings;
use above "Genome"; 
use Test::More tests => 11;

my $e = Genome::Model::Event->get(90005615);
ok(scalar($e), "got a list of events");

for my $io (qw/input output/) {
    my $i;
    if ($io eq 'output') {
        $i = $e->add_output(name => "foooutput", value => "fvaloutput");
    }
    else {
        $i = $e->add_input(name => "fooinput", value => "fvalinput");
    } 
    ok($i, "made an $io");
    ok($i->isa("Genome::Model::Event::" . ucfirst($io)), "$io class is correct: $i");
    is($i->name, "foo$io", "$io name is correct");
    is($i->value, "fval$io", "$io value is correct");
    is($i->event, $e, "$io event is correct");
}
