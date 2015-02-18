#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Test::Exception;

my $pkg = 'Genome::Model::Tools::Joinx::VcfReport';
use_ok($pkg);

my $test_data = File::Spec->catfile($ENV{GENOME_TEST_INPUTS}, "Genome-Model-Tools-Joinx-VcfReport");

my $input_file = File::Spec->catfile($test_data, "input.vcf");
my $expected_site = File::Spec->catfile($test_data, "expected-site.txt");
my $expected_sample = File::Spec->catfile($test_data, "expected-sample.txt");

my $site_out = Genome::Sys->create_temp_file_path;
my $sample_out = Genome::Sys->create_temp_file_path;


my %base_params =(
        per_site_output_file => $site_out,
        per_sample_output_file => $sample_out,
        input_file => $input_file,
        generate_report => 0
        );

my $bad_cmd = $pkg->create(use_version => '1.4', %base_params);
ok($bad_cmd, "created command with invalid version");
dies_ok(sub {$bad_cmd->execute}, "command with invalid version fails to execute");

sub make_versioned_test {
    my $version = shift;

    return sub {
        my $cmd = $pkg->create(
                use_version => $version,
                %base_params,
                );

        ok($cmd, "created_command");
        ok($cmd->execute, "executed command");

        my $site_diff = Genome::Sys->diff_file_vs_file($expected_site, $site_out);
        ok(!$site_diff, "site file is correct") or diag $site_diff;

        my $sample_diff = Genome::Sys->diff_file_vs_file($expected_sample, $sample_out);
        ok(!$sample_diff, "sample file is correct") or diag $sample_diff;
    }
}

for my $minor_version (8..11) {
    my $version = "1.$minor_version";
    subtest "Test version $version" => make_versioned_test($version);
}



done_testing();
