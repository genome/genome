#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;

use_ok('Genome::Model::Tools::Vcf::CreateCrossSampleVcf::BuildClump');

my $get_from_props = Genome::Model::Tools::Vcf::CreateCrossSampleVcf::BuildClump->get(
    backfilled_vcf      => 'a',
    bam_file            => 'b',
    pileup_output_file => 'c',
    sample              => 'd',
    vcf_file            => 'e',
);

my $get_expected_id = '{"backfilled_vcf":"a","bam_file":"b","pileup_output_file":"c","sample":"d","vcf_file":"e"}';
is($get_from_props->id, $get_expected_id, 'id is expected json (get from props)');

my $create_from_props = Genome::Model::Tools::Vcf::CreateCrossSampleVcf::BuildClump->create(
    backfilled_vcf      => 1,
    bam_file            => 2,
    pileup_output_file => 3,
    sample              => 4,
    vcf_file            => 5,
);

my $create_expected_id = '{"backfilled_vcf":"1","bam_file":"2","pileup_output_file":"3","sample":"4","vcf_file":"5"}';
is($create_from_props->id, $create_expected_id, 'id is expected json (create from props)');

my $get_from_id = Genome::Model::Tools::Vcf::CreateCrossSampleVcf::BuildClump->get(
    '{"backfilled_vcf":"alpha","bam_file":"beta","pileup_output_file":"gamma","sample":"delta","vcf_file":"epsilon"}'
);

is($get_from_id->backfilled_vcf      , 'alpha'   , 'backfilled_vcf matches');
is($get_from_id->bam_file            , 'beta'    , 'bam_file matches');
is($get_from_id->pileup_output_file  , 'gamma'   , 'pileup_output_file matches');
is($get_from_id->sample              , 'delta'   , 'sample matches');
is($get_from_id->vcf_file            , 'epsilon' , 'vcf_file matches');

done_testing();
