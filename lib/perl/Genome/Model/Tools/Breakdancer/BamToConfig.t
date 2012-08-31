#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;
use File::Basename;
use File::Compare;
use File::Temp;

my $archos = `uname -a`;
if ($archos !~ /64/) {
    plan skip_all => "Must run from 64-bit machine";
} else {
    plan tests => 8;
}

use_ok( 'Genome::Model::Tools::Breakdancer::BamToConfig');

my $test_input_dir  = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Breakdancer-BamToConfig/';
my $cfg_file    = $test_input_dir . 'breakdancer_config';
my $normal_bam  = $test_input_dir . 'normal.bam';
my $tumor_bam   = $test_input_dir . 'tumor.bam';

my $tmp_dir = File::Temp::tempdir(
    'Genome-Model-Tools-Breakdancer-BamToConfig-XXXXX', 
    DIR     => "$ENV{GENOME_TEST_TEMP}", 
    CLEANUP => 1,
);

my $out_file = $tmp_dir.'/breakdancer_config';

my $bd = Genome::Model::Tools::Breakdancer::BamToConfig->create(
    normal_bam  => $normal_bam,
    tumor_bam   => $tumor_bam,
    output_file => $out_file,
);

$bd->dump_status_messages(1);
ok($bd, 'BamToConfig created ok');
ok($bd->execute(), 'BamToConfig executed ok');

for my $file (glob($tmp_dir."/*insertsize_histogram*"), $out_file) {
    my $base_name = basename $file;
    my $test_file = $test_input_dir . "/$base_name";
    is(compare($file, $test_file), 0, "$base_name output as expected");
}


