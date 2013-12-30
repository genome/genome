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

my $pkg = 'Genome::Model::RnaSeq::Command::RunEdgeR';

use_ok($pkg);
my $num_genes = 4;

my $data_dir = Genome::Utility::Test->data_dir_ok($pkg, "v1");
my $data_dir_n = join("/", $data_dir, "N");
my $data_dir_t = join("/", $data_dir, "T");

my $rnaseq_model_s1_t = Genome::Test::Factory::Model::RnaSeq->setup_object();
my $rnaseq_build_s1_t = Genome::Test::Factory::Build->setup_object(
    model_id        => $rnaseq_model_s1_t->id,
    data_directory  => $data_dir_t,
    status          => "Succeeded",
);
my $rnaseq_model_s1_n = Genome::Test::Factory::Model::RnaSeq->setup_object(
    processing_profile_id  => $rnaseq_model_s1_t->processing_profile_id,
    subject => $rnaseq_model_s1_t->subject,
);
my $rnaseq_build_s1_n = Genome::Test::Factory::Build->setup_object(
    model_id        => $rnaseq_model_s1_n->id,
    data_directory  => $data_dir_n,
    status          => "Succeeded",
);


my $rnaseq_model_s2_t = Genome::Test::Factory::Model::RnaSeq->setup_object(
    processing_profile_id  => $rnaseq_model_s1_t->processing_profile_id,
);
my $rnaseq_build_s2_t = Genome::Test::Factory::Build->setup_object(
    model_id            => $rnaseq_model_s2_t->id,
    data_directory      => $data_dir_t,
    status              => "Succeeded",
);
my $rnaseq_model_s2_n = Genome::Test::Factory::Model::RnaSeq->setup_object(
    processing_profile_id  => $rnaseq_model_s1_t->processing_profile_id,
    subject => $rnaseq_model_s2_t->subject,
);
my $rnaseq_build_s2_n = Genome::Test::Factory::Build->setup_object(
    model_id            => $rnaseq_model_s2_n->id,
    data_directory      => $data_dir_n,
    status              => "Succeeded",
);

my $rnaseq_model_s3_t = Genome::Test::Factory::Model::RnaSeq->setup_object(
    processing_profile_id  => $rnaseq_model_s1_t->processing_profile_id,
);
my $rnaseq_build_s3_t = Genome::Test::Factory::Build->setup_object(
    model_id            => $rnaseq_model_s3_t->id,
    data_directory      => $data_dir_t,
    status              => "Succeeded",
);

subtest "execute" => sub {
    my $output_file = Genome::Sys->create_temp_file_path;

    my $cmd = $pkg->execute(
        normal_models     => [$rnaseq_model_s1_n, $rnaseq_model_s2_n],
        tumor_models      => [$rnaseq_model_s1_t, $rnaseq_model_s2_t],
        counts_per_million => 60000,
        num_samples        => 4,
        output_file        => $output_file,
    );

    my $in = Genome::Sys->open_file_for_reading($output_file);
    my $header = $in->getline;
    chomp $header;
    my @fields = split("\t", $header);
    is_deeply(["", "logFC", "logCPM", "PValue", "FDR", "test.result"], \@fields, "header is as expected");
    my %gene_results;
    my %gene_tests;
    while (my $line = $in->getline) {
        chomp $line;
        my ($name, $x, $y, $pvalue, $adjusted_pvalue, $test) = split("\t", $line);
        $gene_results{$name}{pvalue} = $adjusted_pvalue;
        $gene_results{$name}{test} = $test;
    }

    my @missing = grep {
            my $n = "GENE$_";
            !exists $gene_results{$n}{pvalue} or $gene_results{$n}{pvalue} eq 'NA'
        } 0..$num_genes - 1;

    is(scalar @missing, 0, "Got data for all genes") or diag("Missing: " . Dumper(\@missing));

    ok($gene_results{GENE0}{pvalue} < 0.001,
        "GENE0 is differentially expressed with high confidence");

    is($gene_results{GENE0}{test}, 1,
        "GENE0 passed significance test for DE");

    my @fail = grep { $gene_results{"GENE$_"}{pvalue} < 0.5 } 1..$num_genes - 1;
    is(scalar @fail, 0, "All other genes are p >= 0.5") or diag("Failed: " . Dumper(\@fail));

    @fail = grep { $gene_results{"GENE$_"}{test} != 0 } 1..$num_genes - 1;
    is(scalar @fail, 0, "All other genes fail significance test") or diag("Failed: " . Dumper(\@fail));
};

subtest "execute with unmatched" => sub {
    my $output_file = Genome::Sys->create_temp_file_path;

    my $cmd = $pkg->execute(
        normal_models     => [$rnaseq_model_s1_n, $rnaseq_model_s2_n],
        tumor_models      => [$rnaseq_model_s3_t],
        counts_per_million => 60000,
        num_samples        => 3,
        output_file        => $output_file,
    );

    my $in = Genome::Sys->open_file_for_reading($output_file);
    my $header = $in->getline;
    chomp $header;
    my @fields = split("\t", $header);
    is_deeply(["", "logFC", "logCPM", "PValue", "FDR", "test.result"], \@fields, "header is as expected");
    my %gene_results;
    my %gene_tests;
    while (my $line = $in->getline) {
        chomp $line;
        my ($name, $x, $y, $pvalue, $adjusted_pvalue, $test) = split("\t", $line);
        $gene_results{$name}{pvalue} = $adjusted_pvalue;
        $gene_results{$name}{test} = $test;
    }

    my @missing = grep {
            my $n = "GENE$_";
            !exists $gene_results{$n}{pvalue} or $gene_results{$n}{pvalue} eq 'NA'
        } 0..$num_genes - 1;

    is(scalar @missing, 0, "Got data for all genes") or diag("Missing: " . Dumper(\@missing));

    ok($gene_results{GENE0}{pvalue} < 0.001,
        "GENE0 is differentially expressed with high confidence");

    is($gene_results{GENE0}{test}, 1,
        "GENE0 passed significance test for DE");

    my @fail = grep { $gene_results{"GENE$_"}{pvalue} < 0.5 } 1..$num_genes - 1;
    is(scalar @fail, 0, "All other genes are p >= 0.5") or diag("Failed: " . Dumper(\@fail));

    @fail = grep { $gene_results{"GENE$_"}{test} != 0 } 1..$num_genes - 1;
    is(scalar @fail, 0, "All other genes fail significance test") or diag("Failed: " . Dumper(\@fail));
};

done_testing();
