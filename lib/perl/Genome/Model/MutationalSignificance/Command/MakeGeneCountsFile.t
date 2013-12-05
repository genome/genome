#!/usr/bin/env perl

use above 'Genome';
use Data::Dumper;
use IO::File;
use Test::More;
use Genome::Test::Factory::Model::RnaSeq;
use Genome::Test::Factory::Build;
use Genome::Utility::Test;

use strict;
use warnings;

my $pkg = 'Genome::Model::MutationalSignificance::Command::MakeGeneCountsFile';

use_ok($pkg);
my $data_dir = Genome::Utility::Test->data_dir_ok($pkg, "v1");

my $rna_seq_model = Genome::Test::Factory::Model::RnaSeq->setup_object();
my $rna_seq_build = Genome::Test::Factory::Build->setup_object(model_id => $rna_seq_model->id, data_directory => $data_dir);
my $build_data_directory = $rna_seq_build->data_directory . "/results/digital_expression_result/gene-counts.tsv";
my $build_source = $rna_seq_build->subject->source;

$DB::single=1;
subtest "ok retrieve build information" => sub {
    my @builds = ($rna_seq_build, $rna_seq_build);

    my $build_information = $pkg->_retrieve_build_information("normal", @builds);

    my %expected_information = (
        input_gene_count_files => [$build_data_directory, $build_data_directory],
        subjects => [$build_source, $build_source],
        groups => ["normal", "normal"],
    );

    is_deeply($build_information, \%expected_information, "Build information matches expected build information");
};

subtest "ok create file join command" => sub {
    my @files = ("file1", "file2", "file3", "file4");

    my $cmd = $pkg->_create_file_join_command("output_file", @files);

    is($cmd, "join file1 file2 | join - file3 | join - file4 > output_file", "Join command as expected");
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

subtest "only one group" => sub {
    my $gene_counts_file_path = Genome::Sys->create_temp_file_path;

    my $obj = $pkg->create(
        tumor_rnaseq_builds => [$rna_seq_build, $rna_seq_build],
        normal_rnaseq_builds => [$rna_seq_build, $rna_seq_build],
        gene_counts_file => $gene_counts_file_path,
    );

    my $rv = $obj->execute;

    is($obj->groups_list, "tumor,tumor,normal,normal", "Groups list as expected");
    is($obj->subjects_list, "$build_source,$build_source,$build_source,$build_source", "Subjects list as expected");
};


done_testing();
