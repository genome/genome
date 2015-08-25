#!/usr/bin/env genome-perl
use strict;
use warnings;
BEGIN { $ENV{APP_DBI_NO_COMMIT} = 1 };
use GSCApp;
use above "Genome";
use Test::More tests => 2;

ok($ENV{UR_DBI_NO_COMMIT}, "no-commit transfers from App to UR");

my $o = eval { GSC::Clone->get(10001) };
ok($o, "loading objects works if Genome is loaded after GSCApp")
    || diag("eval error was: $@");
