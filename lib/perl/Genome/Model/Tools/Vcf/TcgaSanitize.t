#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Genome::Utility::Vcf ('diff_vcf_file_vs_file');
use Test::More tests => 13;

use_ok('Genome::Model::Tools::Vcf::TcgaSanitize');

my $test_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Vcf-TcgaSanitize";
my $input_dir = "$test_dir/input.v3";
my $expected_base = "expected.v3";
my $expected_dir = "$test_dir/$expected_base";
my @expected_vcfs = ("$expected_dir/snvs.vcf.gz", "$expected_dir/nontcga.vcf.gz");
my @input_vcfs = ("$input_dir/snvs.vcf.gz", "$input_dir/nontcga.vcf.gz");

for my $index (0..1) {
    my $input_vcf = $input_vcfs[$index];
    my $expected_file = $expected_vcfs[$index];
    my $output_file = Genome::Sys->create_temp_file_path;

    my $command = Genome::Model::Tools::Vcf::TcgaSanitize->create( input_file => $input_vcf, 
        output_file  => $output_file,
        package_for_tcga => 0,
    );

    ok($command, 'Command created');
    my $rv = $command->execute;
    ok($rv, 'Command completed successfully');
    ok(-s $output_file, "output file created");

    my $diff = diff_vcf_file_vs_file($output_file, $expected_file);
    ok(!$diff, "output $output_file matched expected result $expected_file")
        or diag("diff results:\n" . $diff);
}

# Try things with tcga packaging on
my $output_file3 = Genome::Sys->create_temp_directory . "/genome.wustl.edu_CESC.IlluminaGA_DNASeq.Level_2.2.9.0.tar.gz";
my $command2 = Genome::Model::Tools::Vcf::TcgaSanitize->create( input_file => $input_vcfs[0], 
                                                               output_file  => $output_file3,
                                                               package_for_tcga => 1,
                                                             );
ok($command2, 'Command created');
my $rv2 = $command2->execute;
ok($rv2, 'Command completed successfully');
ok(-s $output_file3, "Output zip exists");
ok(-s $output_file3.".md5", "Output zip md5 exists");
