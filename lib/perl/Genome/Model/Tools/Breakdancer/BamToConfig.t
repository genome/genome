#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;
use File::Basename;
use File::Compare qw(compare);
use File::Temp qw();

my $archos = `uname -a`;
if ($archos !~ /64/) {
    plan skip_all => "Must run from 64-bit machine";
} else {
    plan tests => 6;
}

use_ok( 'Genome::Model::Tools::Breakdancer::BamToConfig');

my $test_input_dir  = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Breakdancer-BamToConfig/';
my $cfg_file    = $test_input_dir . 'breakdancer_config';
my $normal_bam  = $test_input_dir . 'normal.bam';
my $tumor_bam   = $test_input_dir . 'tumor.bam';

my $tmp_dir = File::Temp::tempdir(
    'Genome-Model-Tools-Breakdancer-BamToConfig-XXXXX', 
    TMPDIR => 1,
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

for my $file (glob($tmp_dir."/*insertsize_histogram")) {
    my $base_name = basename $file;
    my $test_file = $test_input_dir . "/$base_name";
    is(compare($file, $test_file), 0, "$base_name output as expected");
}

my $orig_config_file = join('/', $test_input_dir, basename $out_file);
my $config_diff = sub {
    my ($line1, $line2) = @_;
    $line1 =~ s/\tmap:[^\t]+\t//;
    $line2 =~ s/\tmap:[^\t]+\t//;
    return $line1 ne $line2;
};
is(compare($orig_config_file, $out_file, $config_diff), 0, "config diffed as expected");
