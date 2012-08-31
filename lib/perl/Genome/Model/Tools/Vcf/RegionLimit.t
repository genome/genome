#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 5;

use_ok('Genome::Model::Tools::Vcf::RegionLimit');

my $test_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Vcf-RegionLimit";
my $expected_base = "expected";
my $expected_dir = "$test_dir/$expected_base";
my $expected_file = "$expected_dir/expected.vcf.gz";

my $output_file = Genome::Sys->create_temp_file_path;
my $input_vcf= "$test_dir/input/snvs.build115632524.vcf.gz";
my $input_bed_regions = "$test_dir/input/feature_list_2.bed";

my $cmd = Genome::Model::Tools::Vcf::RegionLimit->create(   vcf_file => $input_vcf, 
                                                            output_file  => $output_file,
                                                            region_bed_file => $input_bed_regions,
);

ok($cmd, 'Command created');
my $rv = $cmd->execute;
ok($rv, 'Command completed successfully');
ok(-s $output_file, "output file created");

# The files will have a timestamp that will differ. Ignore this but check the rest.
my $expected = `zcat $expected_file |  grep -v fileDate | grep -v region`;
my $output = `zcat $output_file | grep -v fileDate`;
#my $diff = Genome::Sys->diff_text_vs_text($output_file, $expected_file);
my $diff = ($expected eq $output);


ok($diff, 'output matched expected result')
    or diag("diff results:\n" . $diff);
