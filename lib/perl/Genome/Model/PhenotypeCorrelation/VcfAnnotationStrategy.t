#!/usr/bin/env perl

use above 'Genome';
use Test::More;
use Data::Dumper;

use strict;
use warnings;

my $pkg = 'Genome::Model::PhenotypeCorrelation::VcfAnnotationStrategy';

use_ok($pkg);

# List of input examples along with the expected output they should produce
my @tests = (
    # Build specified by source_name/version, specific info fields
    {
        input => 'source dbsnp version 137 info a=b,per-alt:c:d',
        expected => [
            {
                source_name => 'dbsnp',
                source_version => 137,
                info => 'a=b,per-alt:c:d',
            }
        ]
    },
    # Build specified by source_name/version, info fields unspecified
    {
        input => 'source dbsnp version 137',
        expected => [
            {
                source_name => 'dbsnp',
                source_version => 137,
                info => undef,
            }
        ]
    },
    # Build specified by id, specific info fields
    {
        input => 'build 1234 info a=b,per-alt:c:d',
        expected => [ { info => 'a=b,per-alt:c:d', source_build => 1234 } ]
    },
    # Build specified by id, info fields unspecified
    {
        input => 'build 1234',
        expected => [ { info => undef, source_build => 1234 } ]
    },
    # Compound case
    {
        input => join(' ; ', 
            'source dbsnp version 137 info a=b,per-alt:c:d',
            'build 1234',
            'source 1kg-wgs version 20 info x'
            ),
        expected => [
            {
                source_name => 'dbsnp',
                source_version => 137,
                info => 'a=b,per-alt:c:d',
            },
            {
                source_build => 1234,
                info => undef,
            },
            {
                source_name => '1kg-wgs',
                source_version => 20,
                info => 'x',
            },
        ]
    },
);

for my $test (@tests) {
    my $input = $test->{input};
    my $expected = $test->{expected};
    my $obj = $pkg->create(strategy => $input);
    ok($obj, "Created object") or diag("For strategy: $input");
    my $tree = $obj->execute;
    ok($tree, "Execute") or diag("For strategy: $input");
    is_deeply($expected, $tree, "Input string parsed as expected")
        or diag("Expected: " . Dumper($expected) . "Got: " . Dumper($tree));
    is_deeply($tree, $obj->tree, "Object output matches return of execute")
        or diag("Expected: " . Dumper($tree) . "Got: " . Dumper($obj->tree));
}

done_testing();
