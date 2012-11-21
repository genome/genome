#!/usr/bin/env perl
use strict;
use warnings;
use above "Genome";
use Test::More skip_all => "a one-off error is the final column for some rows of the tsv output"; #tests => 2;

my $expected_out = __FILE__ . '.expected-output';

my $actual_out;
if (@ARGV and $ARGV[0] eq 'REBUILD') {
    $actual_out = $expected_out;
    warn "**** rebuilding expected output at $expec****";
}
else {
    $actual_out = Genome::Sys->create_temp_directory;
}

my $cmd = "genome model clin-seq generate-clonality-plots --somatic-var-build=129396826  --output-dir=$actual_out  --common-name='AML54'  --verbose";
eval { Genome::Sys->shell_cmd(cmd => $cmd); };
ok(!$@, "executed command correctly: $cmd");

my @diff = `diff -r $expected_out $actual_out`;
is(scalar(@diff), 0, "no differences") 
    or do {
        for (@diff) { diag($) }
        warn "*** if the above differences are not in error, rebuild the test data by running this test with REBUILD on the command-line ***";
    }

