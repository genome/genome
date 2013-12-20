#!/usr/bin/env genome-perl
use strict;
use warnings;

$ENV{UR_DBI_NO_COMMIT} = 1;

use Test::More;
use above "Genome";

my $class = 'Genome::Config::Rule';
use_ok($class);

my $split_char = '->';

my @rules = Genome::Config::Rule->create_from_hash(_test_hash());

is(scalar(@rules), 2, 'expect to get two Rule objects back');

my %expected_values = reverse %{_test_hash()};
for my $rule (@rules) {
    my $methods = $expected_values{$rule->expected_value};
    my @methods = split($split_char, $methods);
    is(scalar(@{$rule->method_chain}), scalar(@methods), 'expect to get an array ref of correct length');
}

done_testing();

sub _test_hash {
    return {
        "meth1${split_char}meth1${split_char}meth3" => 'three methods',
        'meth4' => 'one method'
    };
}
