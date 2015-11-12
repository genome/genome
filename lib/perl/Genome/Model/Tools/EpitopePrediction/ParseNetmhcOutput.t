#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Utility::Test qw(compare_ok);
use YAML::Syck qw(LoadFile);

my $TEST_DATA_VERSION = 3;
my $class = 'Genome::Model::Tools::EpitopePrediction::ParseNetmhcOutput';
use_ok($class);

my $test_dir = __FILE__ . '.d';

subtest "SNVs file" => sub {
    for my $netmhc_version (qw(3.0 3.4)) {
        my $variant_type = 'snvs';
        for my $output_filter (qw(all top)) {
            subtest "NetMHC version $netmhc_version - $output_filter" => sub {
                test_for_output_type(
                    output_filter => $output_filter,
                    netmhc_version => $netmhc_version,
                    variant_type => $variant_type,
                );
            }
        }
    }
};

subtest "Indel file with frameshifts and inframe indels" => sub {
    my $output_filter = 'all';
    my $netmhc_version = '3.4';
    my $variant_type = 'indels';
    for my $output_filter (qw(all top)) {
        subtest "NetMHC version $netmhc_version - $output_filter" => sub {
            test_for_output_type(
                output_filter => $output_filter,
                netmhc_version => $netmhc_version,
                variant_type => $variant_type,
            );
        };
    }
};

my %match_count_test_cases = (
    K10KM => 'ABCDEFGHIJKLNOPQRSTUV',
    K10KMM => 'ABCDEFGHIJKLNOPQRSTUV',
    K10KMMM => 'ABCDEFGHIJKLNOPQRSTUV',
    K10KMMMMM => 'ABCDEFGHIJKLNOPQRSTUV',
    K10KMMMMMM => 'ABCDEFGHIJKLNOPQRSTUV',
    K10KMMMMMMMMMM => 'ABCDEFGHIJKLNOPQRSTUV',
    K10KMMMMMMMMMMMMMMM => 'ABCDEFGHIJKLNOPQRSTUV',
    K0KM => 'ABCDEFGHIJKLNOPQRSTUV',
    K2KM => 'ABCDEFGHIJKLNOPQRSTUV',
    K5KM => 'ABCDEFGHIJKLNOPQRSTUV',
    K21KM => 'ABCDEFGHIJKLNOPQRSTUV',
    K19KM => 'ABCDEFGHIJKLNOPQRSTUV',
    K16KM => 'ABCDEFGHIJKLNOPQRSTUV',
    KL10K => 'ABCDEFGHIJKLNOPQRSTUVW',
    AB0A => 'ABCDEFGHIJKLNOPQRSTUVW',
    CD2C => 'ABCDEFGHIJKLNOPQRSTUVW',
    FG5F => 'ABCDEFGHIJKLNOPQRSTUVW',
    VW20V => 'ABCDEFGHIJKLNOPQRSTUVW',
    TV19T => 'ABCDEFGHIJKLNOPQRSTUVW',
    RS16R => 'ABCDEFGHIJKLNOPQRSTUVW',
);

subtest "Best Matches" => sub {
    my $best_matches_test_dir = File::Spec->join($test_dir, 'best_matches');

    while (my ($test_case, $wt_sequence) = each (%match_count_test_cases)) {
        subtest $test_case => sub {
            my $wt_file = File::Spec->join($best_matches_test_dir, "$wt_sequence.WT.yaml");
            my $wt_position_data = LoadFile($wt_file);

            my $mt_file = File::Spec->join($best_matches_test_dir, "$test_case.MT.yaml");
            my $mt_match_data = LoadFile($mt_file);

            for my $mt_sequence (keys %{$mt_match_data}) {
                my $expected_data = $mt_match_data->{$mt_sequence};
                my $best_matches = Genome::Model::Tools::EpitopePrediction::ParseNetmhcOutput::best_matches($mt_sequence, $wt_position_data);
                is_deeply($best_matches, $expected_data, "Best matches determined correctly for MT sequences ($mt_sequence)");
            }
        };
    }
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
