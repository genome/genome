#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use strict;
use warnings;

use Test::More tests => 1;

use Genome;

subtest 'regression test for RT 105743' => sub {
    # This is a regression test.  Previously some necessary observers were
    # destroyed during UR::Context->rollback().  "Random" hash ordering would
    # sometimes cause the destruction of observers before other, dependent
    # rollback actions ran.  In this test we rollback twice to not depend on
    # the hash ordering.
    plan tests => 4;
    ok(UR::Context->rollback, 'rolled back to induce destruction of necessary observers');
    my $a_iter = Genome::Disk::Allocation->create_iterator();
    my $a = $a_iter->next;
    ok($a, 'got an allocation');
    ok(UR::Context->rollback, 'rolled back again to induce arhivable side-effects');
    ok(UR::Context->commit, 'commit to demonstrate invalid data for save');
};
