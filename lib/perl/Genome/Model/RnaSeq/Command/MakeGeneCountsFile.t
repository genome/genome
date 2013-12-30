#!/usr/bin/env perl

use above 'Genome';
use Data::Dumper;
use IO::File;
use Test::More;
use Genome::Test::Factory::Model::RnaSeq;
use Genome::Test::Factory::Build;
use Genome::Utility::Test qw/compare_ok/;

use strict;
use warnings;

my $pkg = 'Genome::Model::RnaSeq::Command::MakeGeneCountsFile';

use_ok($pkg);
my $data_dir = Genome::Utility::Test->data_dir_ok($pkg, "v1");

my $rnaseq_model = Genome::Test::Factory::Model::RnaSeq->setup_object();
my $rnaseq_build = Genome::Test::Factory::Build->setup_object(
    model_id        => $rnaseq_model->id,
    data_directory  => $data_dir,
    status          => "Succeeded",
);

my $build_data_directory = $rnaseq_build->data_directory . "/results/digital_expression_result/gene-counts.tsv";
my $build_source = $rnaseq_build->subject->source->common_name;

subtest "write header to file" => sub {
    my @headers = qw(1 2 3 4);

    my $out = Genome::Sys->create_temp_file_path;
    my $rv = $pkg->_write_header_to_file($out, @headers);
    ok($rv, "Subroutine returned 1");
    is(`cat $out`, "GENE\t1\t2\t3\t4\n", "File created as expected");
};

subtest "ok retrieve build information" => sub {
    my @builds = ($rnaseq_build, $rnaseq_build);

    my $build_information = $pkg->_retrieve_build_information("normal", @builds);

    my %expected_information = (
        input_counts_files      => [$build_data_directory, $build_data_directory],
        subjects                => [$build_source, $build_source],
        groups                  => ["normal", "normal"],
        headers                 => ["$build_source-normal", "$build_source-normal"],
    );

    is_deeply($build_information, \%expected_information, "Build information matches expected build information");
};

subtest "ok create file join command" => sub {
    my @files = ("file1", "file2", "file3", "file4");

    my $cmd = $pkg->_create_file_join_command("output_file", @files);

    is($cmd, "join -t '\t' file1 file2 | join -t '\t' - file3 | join -t '\t' - file4 >> output_file", "Join command as expected");
};

subtest "first column on input file not identical" => sub {
    my @files = (
        $data_dir . "/results/digital_expression_result/gene-counts_joined.tsv",
        $data_dir . "/results/digital_expression_result/gene-counts_joined_different.tsv",
    );

    my $rv = 0;
    eval {
        $rv = $pkg->_check_gene_column_identical(@files);
    };
    ok($@, "Different gene columns is an error");
    ok(!$rv);

};

subtest "only one file to join" => sub {
    my @files = ("file1");

    my $rv = 0;
    eval {
        $rv = $pkg->_create_file_join_command(@files);
    };
    ok($@, "Just one file to join is an error");
    ok(!$rv);
};

subtest "execute" => sub {
    my $counts_file_path = Genome::Sys->create_temp_file_path;

    my $obj = $pkg->create(
        normal_models  => [$rnaseq_model, $rnaseq_model],
        tumor_models => [$rnaseq_model, $rnaseq_model],
        counts_file     => $counts_file_path,
    );

    my $rv = $obj->execute;

    is($obj->groups, "tumor,tumor,normal,normal", "Groups list as expected");
    is($obj->subjects, "$build_source,$build_source,$build_source,$build_source", "Subjects list as expected");

    compare_ok(
        $counts_file_path,
        $rnaseq_build->data_directory . "/results/digital_expression_result/gene-counts_joined.tsv",
        "Joined gene counts file as expected"
    );
};


done_testing();
