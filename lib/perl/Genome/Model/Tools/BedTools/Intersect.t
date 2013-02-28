#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 11;

use above 'Genome';

use_ok('Genome::Model::Tools::BedTools');
use_ok('Genome::Model::Tools::BedTools::Intersect');

my $tmp_dir = File::Temp::tempdir('BedTools-Intersect-'.Genome::Sys->username.'-XXXX',CLEANUP => 1, TMPDIR => 1);

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-BedTools-Intersect';

my $bed_file_a = $data_dir .'/a.bed';
my $bed_file_b = $data_dir .'/b.bed';

my $additional_options = [
    { },
    { intersection_type => 'overlaps' },
    { intersection_type => 'a-only' },
];

for my $test_number (0..2) {
    my $expected_file = $data_dir.'/expected/expected.' . $test_number . '.bed';

    my $intersect = Genome::Model::Tools::BedTools::Intersect->create(
        output_file => $tmp_dir .'/test.' . $test_number . '.bed',
        input_file_a => $bed_file_a,
        input_file_b => $bed_file_b,
        input_file_a_format => 'bed',
        %{ $additional_options->[$test_number] },
    );
    
    isa_ok($intersect,'Genome::Model::Tools::BedTools::Intersect');
    ok($intersect->execute,'execute intersectBed command '. $intersect->command_name);

    my $diff = Genome::Sys->diff_file_vs_file($expected_file, $intersect->output_file);
    ok(!$diff, 'expected BED file '. $expected_file .' is identical to '. $intersect->output_file)
        or diag("  diff results:\n" . $diff);
}
