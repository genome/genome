#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Utility::Test qw(compare_ok);
use File::Spec qw();

my $pkg = 'Genome::Model::Tools::Picard::CleanSam';

use_ok($pkg);

my $data_dir = File::Spec->catfile(Genome::Config::get('test_inputs'), 'Genome-Model-Tools-Picard-CleanSam');
my $sam_file = File::Spec->catfile($data_dir, 'dirty.sam');
my $expected_file = File::Spec->catfile($data_dir, 'expected.sam');
my $output_file = Genome::Sys->create_temp_file_path . ".sam";

my $obj = $pkg->create(
    input_file => $sam_file,
    output_file => $output_file,
    use_version => "1.123",
    );

ok($obj, 'created command');
ok($obj->execute, 'executed command');

compare_ok($expected_file, $output_file, diag => 1, filters => ['^@.*']);
done_testing();
