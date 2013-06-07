#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 5;

use_ok('Genome::Model::Tools::Vcf::VcfMergeAnnotation');

my $test_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Vcf-VcfMergeAnnotation";
my $output_file = Genome::Sys->create_temp_file_path;
my $input_dir = "$test_dir";
my $input_vcf = "$input_dir/annotation.vcf";
my $input_vep = "$input_dir/vep.output";
my $expected_file = "$input_dir/AnnotationMerged.Expected.vcf";

my $command= Genome::Model::Tools::Vcf::VcfMergeAnnotation->create(
    output_file_name => $output_file,
    vcf_file_name => $input_vcf,
    vep_file_name => $input_vep,
);

ok($command, 'Command created');
my $rv = $command->execute;
ok($rv, 'Command completed successfully');
ok(-s $output_file, "output file created");

# The files will have a timestamp that will differ. Ignore this but check the rest.
my $diff = Genome::Sys->diff_file_vs_file($output_file, $expected_file);
ok($diff, 'output matched expected result')
