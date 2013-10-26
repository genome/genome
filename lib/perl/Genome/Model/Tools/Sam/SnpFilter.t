#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use Test::More tests => 4;
use File::Temp;
use File::Copy;
use File::Compare;

BEGIN {
    use_ok('Genome::Model::Tools::Sam::SnpFilter');
}

my $root_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sam/SnpFilter';

my $tmp_dir  = File::Temp::tempdir(
    "SnpFilter_XXXXXX", 
    TMPDIR => 1,
    CLEANUP => 1,
);

my $snp_file   = "$root_dir/test.snp";
my $out_file   = "$tmp_dir/test.snp.sam_SNPfilter";
my $indel_file = "$root_dir/test.indel";
my $ori_file   = "$root_dir/test.snp.filtered.ori";

my $filter = Genome::Model::Tools::Sam::SnpFilter->create(
    snp_file   => $snp_file,                                                      
    out_file   => $out_file,
    indel_file => $indel_file,
);

isa_ok($filter,'Genome::Model::Tools::Sam::SnpFilter');
ok($filter->execute,'executed ok');

cmp_ok(compare($out_file, $ori_file), '==', 0, 'Sam SNPfilter file was created ok');
