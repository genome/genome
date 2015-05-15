#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 8;

use above 'Genome';
use Genome::Utility::Test qw(compare_ok);

my $pkg = 'Genome::Model::Tools::BedTools::PairToPair';
use_ok($pkg);

my $version = "v1";
my $data_dir = Genome::Utility::Test->data_dir_ok($pkg, $version) or die;

my $bedpe_file_a = $data_dir .'/a.bedpe';
my $bedpe_file_b = $data_dir .'/b.bedpe';

my $additional_options = [
    { },
    { intersection_type => 'notboth' },
];

for my $test_number (0..1) {
    my $expected_file = $data_dir.'/expected.' . $test_number . '.bedpe';

    my $output_file = Genome::Sys->create_temp_file_path;
    my $command = Genome::Model::Tools::BedTools::PairToPair->create(
        output_file => $output_file,
        input_file_a => $bedpe_file_a,
        input_file_b => $bedpe_file_b,
        slop => 0,
        ignore_strand => 1,
        use_version => '2.17.0',
        %{ $additional_options->[$test_number] },
    );

    isa_ok($command,'Genome::Model::Tools::BedTools::PairToPair');
    ok($command->execute,'execute command '. $command->command_name);

    compare_ok($expected_file, $output_file, "bedpe #$test_number is correct");
}
