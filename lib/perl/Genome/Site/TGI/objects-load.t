#!/usr/bin/env genome-perl
use strict;
use warnings;
BEGIN { $ENV{UR_DBI_NO_COMMIT} = 1 };

use above "Genome";
use Test::More tests => 4;

ok(!$GSCApp::{"BEGIN"}, "the GSC namespace is not imported with Genome");

my $v1 = \&UNIVERSAL::isa;

my $o = GSC::Clone->get(10001);
ok($o, "got an object from the GSC namespace");

my $v2 = \&UNIVERSAL::isa;

is($v2,$v1,"the UNIVERSAL::isa method has the same value after using GSCApp");

is($ENV{APP_DBI_NO_COMMIT},1, "no-commit transfers UR-to-App");

