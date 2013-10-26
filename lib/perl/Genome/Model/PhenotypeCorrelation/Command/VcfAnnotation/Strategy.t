#!/usr/bin/env perl

use above 'Genome';
use Test::More;
use Data::Dumper;
use Text::Balanced qw/extract_variable/;

use strict;
use warnings;

my $pkg = 'Genome::Model::PhenotypeCorrelation::Command::VcfAnnotation::Strategy';

use_ok($pkg);

# List of input examples along with the expected output they should produce
my @tests = (
    # Build specified by source_name/version, specific info fields
    {
        input => 'joinx {source_name: dbsnp, source_version: 137, info_fields: "a=b,per-alt:c:d"}',
        expected => [
            {
                type => "joinx",
                params => [{
                    source_name => 'dbsnp',
                    source_version => 137,
                    info_fields => 'a=b,per-alt:c:d',
                }],
            }
        ]
    },
    # Build specified by source_name/version, info fields unspecified
    {
        input => 'joinx {source_name: dbsnp, source_version: 137}',
        expected => [
            {
                type => "joinx",
                params => [{
                    source_name => 'dbsnp',
                    source_version => 137,
                }],
            }
        ]
    },
    ## Build specified by id, specific info fields
    {
        input => 'joinx {source_build: 1234, info_fields: "a=b,per-alt:c:d"}',
        expected => [
            {
                type => "joinx",
                params => [{
                    info_fields => 'a=b,per-alt:c:d',
                    source_build => 1234
                }],
            }
        ]
    },
    ## Build specified by id, info fields unspecified
    {
        input => 'joinx {source_build: 1234}',
        expected => [
            {
                type => "joinx",
                params => [{
                    source_build => 1234
                }],
            }
        ]
    },
    ## Compound case
    {
        input => join(' | ',
            'vep {ensembl_annotation_build: 12}',
            'joinx [{source_name: dbsnp, source_version: 137, info_fields: "a=b,per-alt:c:d"}, {source_build: 1234}, {source_name: 1kg-wgs, source_version: 20, info_fields: x}]',
            ),
        expected => [
            {
                type => "vep",
                params => {
                    ensembl_annotation_build => 12,
                }
            },
            {
                type => "joinx",
                params => [{
                    source_name => 'dbsnp',
                    source_version => 137,
                    info_fields => 'a=b,per-alt:c:d',
                },
                {
                    source_build => 1234,
                },
                {
                    source_name => '1kg-wgs',
                    source_version => 20,
                    info_fields => 'x',
                }]
            }
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
        or diag("Input: $input\nExpected: " . Dumper($expected) . "Got: " . Dumper($tree));
    is_deeply($tree, $obj->tree, "Object output matches return of execute")
        or diag("Input: $input\nExpected: " . Dumper($tree) . "Got: " . Dumper($obj->tree));
}

done_testing();
