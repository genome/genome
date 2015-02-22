#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use Test::More tests => 2;

use_ok('Genome::Config::Tag');

my $tag = Genome::Config::Tag->create(
    name => 'testing_tag_creation',
    description => 'A phony tag just to test instantiation of tags'
);
isa_ok($tag, 'Genome::Config::Tag');

