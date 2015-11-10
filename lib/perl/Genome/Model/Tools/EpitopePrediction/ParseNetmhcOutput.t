#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Utility::Test qw(compare_ok);

my $TEST_DATA_VERSION = 3;
my $class = 'Genome::Model::Tools::EpitopePrediction::ParseNetmhcOutput';
use_ok($class);

my $test_dir = __FILE__ . '.d';

for my $netmhc_version (qw(3.0 3.4)) {
    my $variant_type = 'snvs';
    for my $output_filter (qw(all top)) {
        subtest "NetMHC version $netmhc_version - $variant_type - $output_filter" => sub {
            test_for_output_type(
                output_filter => $output_filter,
                netmhc_version => $netmhc_version,
                variant_type => $variant_type,
            );
        }
    }
}

subtest "Indel file with frameshifts and inframe indels" => sub {
    my $output_filter = 'all';
    my $netmhc_version = '3.4';
    my $variant_type = 'indels';
    test_for_output_type(
        output_filter => $output_filter,
        netmhc_version => $netmhc_version,
        variant_type => $variant_type,
    );
};

sub test_for_output_type {
    my %params = @_;

    my $netmhc_version = $params{netmhc_version};
    my $variant_type = $params{variant_type};
    my $output_filter = $params{output_filter};

    my $expected_output = File::Spec->join($test_dir, "parsed_file.$netmhc_version.$variant_type.$output_filter");
    delete $params{variant_type};

    my $key_file = File::Spec->join($test_dir, "key_file.$variant_type.txt");
    my $netmhc_file = File::Spec->join($test_dir, "netmhc_file.$netmhc_version.$variant_type.xls");
    my $output_dir = Genome::Sys->create_temp_directory;

    my $cmd = $class->create(
        %params,
        netmhc_file => $netmhc_file,
        key_file => $key_file,
        output_directory => $output_dir,
    );
    ok($cmd, "Command created successfully");

    ok($cmd->execute, "Command executes successfully");

    compare_ok($cmd->parsed_file, $expected_output, "Output file is as expected");
}

done_testing();
