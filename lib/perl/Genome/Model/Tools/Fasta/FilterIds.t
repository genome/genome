#!/usr/bin/env genome-perl
use strict;
use warnings;
use above 'Genome';
use Test::More tests => 4;
use_ok("Genome::Model::Tools::Fasta::FilterIds");
my $in = __FILE__ . ".in.fa";
my $out = Genome::Sys->create_temp_file_path("FilterIds.t.out.fa");
my $expected = __FILE__ . ".expected.fa";

if ($ARGV[0] eq 'REBUILD') {
    warn "rebuilding expected output at $expected";
    $out = $expected;
}

eval {
    Genome::Model::Tools::Fasta::FilterIds->execute(
        input_filename => $in,
        output_filename => $out,
        whitelist_regex => "[A-Z]",
        blacklist_regex => "B",
        verbose => 1,
    );
};

ok(!$@, "no exceptions during execution");
ok(-e $out, "output file exists");

my @diff = `diff $out $expected`;
is(scalar(@diff), 0, "no differences from expected");


