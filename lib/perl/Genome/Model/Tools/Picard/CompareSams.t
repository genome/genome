#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use Test::More tests => 9;

use_ok('Genome::Model::Tools::Picard::CompareSams');

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Picard-CompareSams/';

my $file1 = $data_dir . 'normal.tiny.bam';
my $file2 = $data_dir . 'tumor.tiny.bam';

#Test default version
my $output_1 = Genome::Sys->create_temp_file_path();

my $command_1 = Genome::Model::Tools::Picard::CompareSams->create(
    input_file_1 => $file1,
    input_file_2 => $file2,
    output_file => $output_1,
);

ok($command_1, 'created command for default picard version');
ok($command_1->execute(), 'executed command');
ok(-s $output_1, 'produced some output');

#Test an old version
my $output_2 = Genome::Sys->create_temp_file_path();

my $command_2 = Genome::Model::Tools::Picard::CompareSams->create(
    input_file_1 => $file1,
    input_file_2 => $file2,
    output_file => $output_2,
    use_version => 'r116',
);

ok($command_2, 'created command for picard version r116');
ok(! $command_2->execute(), 'could not execute command for unsupported version');

#Test a newer version
my $output_3 = Genome::Sys->create_temp_file_path();

my $command_3 = Genome::Model::Tools::Picard::CompareSams->create(
    input_file_1 => $file1,
    input_file_2 => $file2,
    output_file => $output_3,
    use_version => '1.25', #test new classpath location
);

ok($command_3, 'created command for picard version 1.25');
ok($command_3->execute(), 'executed command');
ok(-s $output_3, 'produced some output');
