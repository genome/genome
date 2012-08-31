#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use File::Compare;
use Test::More tests => 4;

use_ok('Genome::Model::Tools::BedTools::Slop');

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-BedTools-Slop';

my $test_bed = $data_dir .'/test.bed';
my $genome_file = $data_dir .'/all_sequences.genome';
my $expected_bed = $data_dir .'/expected_250bp_slop.bed';

my $tmp_dir = File::Temp::tempdir('BedTools-Slop-'.Genome::Sys->username.'-XXXX',DIR => "$ENV{GENOME_TEST_TEMP}",CLEANUP => 1);
my $output_bed = $tmp_dir .'/slop.bed';

my $slop = Genome::Model::Tools::BedTools::Slop->create(
    use_version => '2.14.3',
    input_file => $test_bed,
    genome_file => $genome_file,
    output_file => $output_bed,
    both => '250',
);
isa_ok($slop,'Genome::Model::Tools::BedTools::Slop');
ok($slop->execute,'execute command '. $slop->command_name);

ok( (compare($slop->output_file,$expected_bed) == 0),'slop BED are identical');

exit;
