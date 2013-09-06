#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use Test::More tests => 4;
use File::Temp;
use File::Compare;

BEGIN {
    use_ok('Genome::Model::Tools::Sam::SnpSanitizer');
}

my $root_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sam/SnpSanitizer';

my $tmp_dir  = File::Temp::tempdir(
    "SnpSanitizer_XXXXXX", 
    TMPDIR => 1,
    CLEANUP => 1,
);

my $snp_file = "$root_dir/test_snp.out";
my $out_file = "$tmp_dir/test_snp.out.sanitize";
my $ori_file = "$root_dir/test_snp.out.sanitize";

my $sani = Genome::Model::Tools::Sam::SnpSanitizer->create(
    snp_file => $snp_file,                                                      
    out_file => $out_file,
);

isa_ok($sani,'Genome::Model::Tools::Sam::SnpSanitizer');
ok($sani->execute,'executed ok');

cmp_ok(compare($out_file, $ori_file), '==', 0, 'Sam SnpSanitizer file was created ok');
