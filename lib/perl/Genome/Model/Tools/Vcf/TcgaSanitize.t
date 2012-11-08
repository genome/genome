#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Genome::Utility::Vcf ('diff_vcf_file_vs_file');
use Test::More tests => 5;

use_ok('Genome::Model::Tools::Vcf::Backfill');

my $test_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Vcf-TcgaSanitize";
my $input_dir = "$test_dir/input.v1";
my $expected_base = "expected.v1";
my $expected_dir = "$test_dir/$expected_base";
my $expected_file = "$expected_dir/snvs.vcf.gz";

my $output_file = Genome::Sys->create_temp_file_path;
my $input_vcf= "$input_dir/snvs.vcf.gz";

my $command = Genome::Model::Tools::Vcf::TcgaSanitize->create( input_file => $input_vcf, 
                                                               output_file  => $output_file,
                                                               package_for_tcga => 0,
                                                             );

ok($command, 'Command created');
my $rv = $command->execute;
ok($rv, 'Command completed successfully');
ok(-s $output_file, "output file created");

my $diff = diff_vcf_file_vs_file($output_file, $expected_file);
ok(!$diff, 'output matched expected result')
    or diag("diff results:\n" . $diff);

=cut
my $output_file2 = Genome::Sys->create_temp_directory . "/genome.wustl.edu_CESC.IlluminaGA_DNASeq.Level_2.2.9.0.vcf";
my $command2 = Genome::Model::Tools::Vcf::TcgaSanitize->create( input_file => $input_vcf, 
                                                               output_file  => $output_file2,
                                                               package_for_tcga => 1,
                                                             );
ok($command2, 'Command created');
my $rv2 = $command2->execute;
ok($rv2, 'Command completed successfully');
my $expected_output_zip = $output_file2 =~ s/vcf/tar.gz/;
ok(-s $expected_output_zip, "Output zip exists");
