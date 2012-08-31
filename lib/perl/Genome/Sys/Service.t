#!/usr/bin/env genome-perl
use strict;
use warnings;
use above 'Genome';
use Test::More tests => 3;

use_ok("Genome::Sys::Service");

my @services = Genome::Sys::Service->get();
ok(scalar(@services) > 0, "got " . scalar(@services) . " services");

my @bad = grep { not ($_->isa("Genome::Sys::Service") and $_->isa("UR::Singleton")) } @services;
is(scalar(@bad), 0, "got zero bad services, as expected");


