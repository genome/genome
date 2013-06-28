#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 3;

use_ok('Genome::Model::Build::ImportedReferenceSequence');

SKIP: {
    unless (-x "/gsc/bin/perl") {
        skip "No /gsc/bin/perl available", 2;
    }

    # This test was added due to a production bug.
    # If it turns out to be fragile, you can just delete it.
    system("/gsc/bin/perl", "-e",
        "require GSCApp; use above 'Genome'; require Genome::Model::ReferenceSequence;");
    my $exit_code = $? >> 8;
    my $signal = $? & 127;
    is($exit_code, 0,
        'nonzero exit code loading ImportedReferenceSequence with GSCApp');
    is($signal, 0,
        'nonzero signal loading ImportedReferenceSequence with GSCApp');
};
