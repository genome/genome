#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::Build::GenotypeMicroarray') or die;
use_ok('Genome::Model::GenotypeMicroarray::Test') or die;

my $build = Genome::Model::GenotypeMicroarray::Test->build;
ok($build, 'got genotype microarray build');

# No variation list [dbsnp_build] build
my $variation_list_build = $build->dbsnp_build;
ok(!$build->dbsnp_build(undef), 'unset variation list build on gm build');
my @tags = $build->validate_dbsnp_build;
ok(@tags, 'validate_dbsnp_build returned errors') or die;
is($tags[0]->desc, 'No DB Snp build specified for build!', 'correct error message');
ok($build->dbsnp_build($variation_list_build), '[re]set variation list build on gm build');

# Variation list build [dbsnp_build] does not have VCF
my $output_dir = $variation_list_build->snv_result->output_dir;
is($variation_list_build->snv_result->output_dir('.'), '.', 'set output dir on snv result to "."');
@tags = $build->validate_dbsnp_build;
ok(@tags, 'validate_dbsnp_build returned errors') or die;
like($tags[0]->desc, qr/ does not have a SNVS VCF!$/, 'correct error message');
is($variation_list_build->snv_result->output_dir($output_dir), $output_dir, '[re]set output dir on snv result');

done_testing()
