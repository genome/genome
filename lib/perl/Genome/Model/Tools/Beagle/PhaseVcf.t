#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Genome::Utility::Vcf ('diff_vcf_file_vs_file');
use Test::More tests => 13;

use_ok('Genome::Model::Tools::Beagle::PhaseVcf');

my $test_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Beagle-PhaseVcf";
my $input_dir = "$test_dir/input.v1";
my $expected_base = "expected.v1";
#my $expected_dir = "$test_dir/$expected_base";
#my $expected_vcf = "$expected_dir/expected.vcf.gz";
my $input_vcf = "$input_dir/input.vcf.gz";
my $input_ped = "input_dir/input.ped";
my $output_file = Genome::Sys->create_temp_file_path;

my $command = Genome::Model::Tools::Beagle::PhaseVcf->create( input_file => $input_vcf, 
    output_file  => $output_file,
    ped_file => $input_ped,
);

ok($command, 'Command created');
my $rv = $command->execute;
ok($rv, 'Command completed successfully');
ok(-s $output_file, "output file created");

#    my $diff = diff_vcf_file_vs_file($output_file, $expected_file);
#    ok(!$diff, "output $output_file matched expected result $expected_file")
#        or diag("diff results:\n" . $diff);
}

