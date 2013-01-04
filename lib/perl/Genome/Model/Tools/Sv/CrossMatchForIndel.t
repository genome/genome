#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;
use File::Slurp;

BEGIN {
    my $archos = `uname -a`;
    if ($archos !~ /64/) {
        plan skip_all => "Must run from 64-bit machine";
    } 
    else {
        plan tests => 4;
    }
};

use_ok( 'Genome::Model::Tools::Sv::CrossMatchForIndel');

my $test_input_dir  = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sv-CrossMatchForIndel/';
my $cm_file  = $test_input_dir . 'cm.file';
my $ref_seq  = $test_input_dir . 'ref.seq';
my $ctg_file = $test_input_dir . 'contig.seq';
my $out_file = $test_input_dir . 'out.file';

my $cm = Genome::Model::Tools::Sv::CrossMatchForIndel->create(
    #output_file          => $test_out_file,
    cross_match_file     => $cm_file,
    local_ref_seq_file   => $ref_seq,
    assembly_contig_file => $ctg_file,
    per_sub_rate         => 0.02,
    ref_start_pos        => '16_8646775_16_8646775_INS_93_+-',
);

ok($cm, 'created CrossMatchForIndel object ok');

my $out_str = $cm->execute();
ok($out_str, 'executed CrossMatchForIndel object OK');

my $ori_out_str = read_file($out_file);
cmp_ok($out_str, 'eq', $ori_out_str, 'output string is generated as expected');
