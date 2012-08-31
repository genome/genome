#!/usr/bin/env genome-perl
use strict;
use warnings;
use above "Genome";
use Test::More tests => 6;

ok(!$GSCApp::{"BEGIN"}, "the GSC namespace is not imported with Genome");

my $v1 = \&UNIVERSAL::isa;

eval "use GSCApp";
ok(!$@, "used GSCApp does nothing but does not crash before really using any GSC modules");

my $v2 = \&UNIVERSAL::isa;
is($v2,$v1,"the UNIVERSAL::isa method has the same value after using GSCApp");

my $o = eval { GSC::Clone->get(10001) };
ok($o, "finds objects int the GSC namespace") or diag $@;

my $v3 = \&UNIVERSAL::isa;
is($v3,$v1,"the UNIVERSAL::isa method has the same value after using GSCApp");

eval "use GSCApp";
ok(!$@, "used GSCApp does nothing but does not crash after actually using GSC modules");

