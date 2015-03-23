#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Utility::Test qw(compare_ok);
use File::Spec qw();

my $pkg = 'Genome::Model::Tools::Picard::CollectWgsMetrics';

use_ok($pkg);

my $data_dir = Genome::Utility::Test->data_dir_ok($pkg);

my $bam_file = File::Spec->join($data_dir, 'sorted.bam');
my $ref_file = File::Spec->join($data_dir, 'small.fa');
my $output_file = Genome::Sys->create_temp_file_path;

my $obj = $pkg->create(
    input_file => $bam_file,
    output_file => $output_file,
    reference_sequence => $ref_file,
    use_version => 1.123
);
ok($obj, 'created command');
ok($obj->execute, 'executed command');

my $expected_file = File::Spec->catfile($data_dir, 'expected.txt');
compare_ok($expected_file, $output_file, filters => ['^#.*']);

done_testing();
