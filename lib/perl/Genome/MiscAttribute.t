#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More tests => 7;

my $x = Genome::InstrumentData::Solexa->create(id => -1);
ok($x, "created a solexa lane object");

$x->full_path("/my/path");
is($x->full_path,'/my/path', "set the full path on the main object");

my $a = $x->add_attribute(
    attribute_label => "full_path2",
    attribute_value => "/my/path"
);
ok($a, "added an attribute the 'hard way'");
is($a->attribute_label, "full_path2", "property name is correct");
is($a->attribute_value, "/my/path", "vaue is correct");
is($a->instrument_data,$x, "entity object is correct");
is($a->instrument_data_id, $x->id, "entity id is correct");

#my @d = Genome::InstrumentData::Solexa->get(full_path => "/my/path");
#is(scalar(@d),1, "got expected object returned by query");
