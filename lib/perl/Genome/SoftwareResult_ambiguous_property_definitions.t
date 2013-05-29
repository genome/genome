#!/usr/bin/env perl

use strict;
use warnings;
use Test::More tests => 2;

use above 'Genome';

my $classname = 'BadSoftwareResult';
my $class = eval {
    UR::Object::Type->define(
        class_name => $classname,
        is => 'Genome::SoftwareResult',
        has => [
            ambiguous_property => {
                is => 'Text',
                # If a subclass does not specifiy an alternate table name then UR assumes the
                # properties do not need to be saved to the superclasses data source. So this
                # test is to try to force properties to be set to either is_input, is_param, or
                # is_transient (for non-calculated or non-delegated properties).
                #
                # We had ~300k objects affected by this, see RT #88011.
            },
        ],
    );
};
my $error = $@;
ok(!defined($class), qq(failed to define $classname));
like($error, qr/ambiguous_property/, 'got an error about ambiguous_property');
