#!/usr/bin/env genome-perl

use strict;
use warnings;
use above 'Genome';
use Test::More tests => 2;
use Genome::Utility::Test qw(command_execute_ok);


my $class = 'Genome::Model::Tools::ListReferenceSequences';
use_ok($class);


do {
    use Genome::Model::Build;
    no warnings 'once';
    local *Genome::Model::Build::get = sub {
        my $class = shift;
        my @got = UR::Object::get($class, @_);
        if (@got > 5) {
            note 'Limiting to 5 results for test.';
            @got = @got[0..4];
        }
        return @got;
    };
    use warnings 'once';

    my $cmd = $class->create();
    command_execute_ok($cmd, 'ListReferenceSequences ran OK');
};
