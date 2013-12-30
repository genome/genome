#!/usr/bin/env perl

use above 'Genome';
use Data::Dumper;
use IO::File;
use Test::More;
use Genome::Test::Factory::Model::RnaSeq;
use Genome::Test::Factory::Model::ClinSeq;
use Genome::Test::Factory::Build;
use Genome::Utility::Test qw/compare_ok/;

use strict;
use warnings;

my $pkg = 'Genome::Model::MutationalSignificance::Command::MakeGeneCountsFile';

use_ok($pkg);
my $data_dir = Genome::Utility::Test->data_dir_ok($pkg, "v2");

my $rnaseq_model = Genome::Test::Factory::Model::RnaSeq->setup_object();
my $rnaseq_build = Genome::Test::Factory::Build->setup_object(
    model_id        => $rnaseq_model->id,
    data_directory  => $data_dir,
    status          => "Succeeded",
);

my $clinseq_model = Genome::Test::Factory::Model::ClinSeq->setup_object(
    tumor_rnaseq_model  => $rnaseq_model,
    normal_rnaseq_model => $rnaseq_model,
);

my $build_data_directory = $rnaseq_build->data_directory . "/results/digital_expression_result/gene-counts.tsv";
my $build_source = $rnaseq_build->subject->source->common_name;

subtest "execute" => sub {
    my $counts_file_path = Genome::Sys->create_temp_file_path;

    my $obj = $pkg->create(
        clinseq_models  => [$clinseq_model, $clinseq_model],
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
