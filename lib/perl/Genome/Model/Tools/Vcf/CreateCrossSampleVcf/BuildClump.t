#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;

use_ok('Genome::Model::Tools::Vcf::CreateCrossSampleVcf::BuildClump');

my $get_from_props = Genome::Model::Tools::Vcf::CreateCrossSampleVcf::BuildClump->get(
    backfilled_vcf      => 'a',
    bam_file            => 'b',
    build_id            => 'c',
    filtered_vcf        => 'd',
    pileup_output_file => 'e',
    sample              => 'f',
    vcf_file            => 'g',
);

my $get_expected_id = '{"backfilled_vcf":"a","bam_file":"b","build_id":"c","filtered_vcf":"d","pileup_output_file":"e","sample":"f","vcf_file":"g"}';
is($get_from_props->id, $get_expected_id, 'id is expected json (get from props)');

my $create_from_props = Genome::Model::Tools::Vcf::CreateCrossSampleVcf::BuildClump->create(
    backfilled_vcf      => 1,
    bam_file            => 2,
    build_id            => 3,
    filtered_vcf        => 4,
    pileup_output_file  => 5,
    sample              => 6,
    vcf_file            => 7,
);

my $create_expected_id = '{"backfilled_vcf":"1","bam_file":"2","build_id":"3","filtered_vcf":"4","pileup_output_file":"5","sample":"6","vcf_file":"7"}';
is($create_from_props->id, $create_expected_id, 'id is expected json (create from props)');

my $get_from_id = Genome::Model::Tools::Vcf::CreateCrossSampleVcf::BuildClump->get(
    '{"backfilled_vcf":"alpha","bam_file":"beta","build_id":"gamma","filtered_vcf":"delta","pileup_output_file":"epsilon","sample":"zeta","vcf_file":"eta"}'
);

is($get_from_id->backfilled_vcf      , 'alpha'   , 'backfilled_vcf matches');
is($get_from_id->bam_file            , 'beta'    , 'bam_file matches');
is($get_from_id->build_id            , 'gamma'   , 'bam_file matches');
is($get_from_id->filtered_vcf        , 'delta'   , 'bam_file matches');
is($get_from_id->pileup_output_file  , 'epsilon' , 'pileup_output_file matches');
is($get_from_id->sample              , 'zeta'    , 'sample matches');
is($get_from_id->vcf_file            , 'eta'     , 'vcf_file matches');

done_testing();
