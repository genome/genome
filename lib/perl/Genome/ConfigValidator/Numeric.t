#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 5;

my $class = 'Genome::ConfigValidator::Numeric';
use_ok($class);

my $validator = $class->new();
my @tests = (
    [123, 1],
    ['123', 1],
    ['abc', 0],
    ['123 abc', 0],
);
for my $test (@tests) {
    my ($value, $expected) = @{$test};
    my $message = sprintf(q('%s' should %s check), $value, ($expected ? 'pass' : 'fail'));
    is($validator->check($value), $expected, $message);
}
