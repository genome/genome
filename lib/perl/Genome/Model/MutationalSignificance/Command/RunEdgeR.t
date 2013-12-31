#!/usr/bin/env perl

use above 'Genome';
use Data::Dumper;
use IO::File;
use Test::More;
use Genome::Test::Factory::Model::RnaSeq;
use Genome::Test::Factory::Model::ClinSeq;
use Genome::Test::Factory::Model::ImportedAnnotation;
use Genome::Test::Factory::Build;
use Genome::Utility::Test qw/compare_ok/;

use strict;
use warnings;

my $pkg = 'Genome::Model::MutationalSignificance::Command::RunEdgeR';

use_ok($pkg);
my $num_genes = 4;

my $data_dir_t = Genome::Utility::Test->data_dir_ok($pkg, "T");
my $data_dir_n = Genome::Utility::Test->data_dir_ok($pkg, "N");

my $annot_model = Genome::Test::Factory::Model::ImportedAnnotation->setup_object();
my $annot_build = Genome::Test::Factory::Build->setup_object(
    model_id => $annot_model->id,
    status   => "Succeeded",
);

my $rnaseq_model_s1_t = Genome::Test::Factory::Model::RnaSeq->setup_object();
my $rnaseq_build_s1_t = Genome::Test::Factory::Build->setup_object(
    model_id        => $rnaseq_model_s1_t->id,
    data_directory  => $data_dir_t,
    status          => "Succeeded",
    annotation_build => $annot_build,
);
my $rnaseq_model_s1_n = Genome::Test::Factory::Model::RnaSeq->setup_object(
    processing_profile_id  => $rnaseq_model_s1_t->processing_profile_id,
    subject => $rnaseq_model_s1_t->subject,
);
my $rnaseq_build_s1_n = Genome::Test::Factory::Build->setup_object(
    model_id        => $rnaseq_model_s1_n->id,
    data_directory  => $data_dir_n,
    status          => "Succeeded",
    annotation_build => $annot_build,
);


my $rnaseq_model_s2_t = Genome::Test::Factory::Model::RnaSeq->setup_object(
    processing_profile_id  => $rnaseq_model_s1_t->processing_profile_id,
);
my $rnaseq_build_s2_t = Genome::Test::Factory::Build->setup_object(
    model_id            => $rnaseq_model_s2_t->id,
    data_directory      => $data_dir_t,
    status              => "Succeeded",
    annotation_build => $annot_build,
);
my $rnaseq_model_s2_n = Genome::Test::Factory::Model::RnaSeq->setup_object(
    processing_profile_id  => $rnaseq_model_s1_t->processing_profile_id,
    subject => $rnaseq_model_s2_t->subject,
);
my $rnaseq_build_s2_n = Genome::Test::Factory::Build->setup_object(
    model_id            => $rnaseq_model_s2_n->id,
    data_directory      => $data_dir_n,
    status              => "Succeeded",
    annotation_build => $annot_build,
);

my $clinseq_model_s1 = Genome::Test::Factory::Model::ClinSeq->setup_object(
    tumor_rnaseq_model  => $rnaseq_model_s1_t,
    normal_rnaseq_model => $rnaseq_model_s1_n,
);

my $clinseq_model_s2 = Genome::Test::Factory::Model::ClinSeq->setup_object(
    tumor_rnaseq_model  => $rnaseq_model_s2_t,
    normal_rnaseq_model => $rnaseq_model_s2_n,
);

subtest "execute" => sub {
    my $output_file = Genome::Sys->create_temp_file_path;

    my $cmd = $pkg->execute(
        clinseq_models     => [$clinseq_model_s1, $clinseq_model_s2],
        counts_per_million => 60000,
        num_samples        => 4,
        output_file        => $output_file,
    );

    my $in = Genome::Sys->open_file_for_reading($output_file);
    my $header = $in->getline;
    chomp $header;
    my @fields = split("\t", $header);
    is_deeply(["", "genes", "logFC", "logCPM", "LR", "PValue", "FDR", "test.result"], \@fields, "header is as expected");
    my %gene_results;
    my %gene_tests;
    while (my $line = $in->getline) {
        chomp $line;
        my ($coeff, $name, $x, $y, $z, $pvalue, $adjusted_pvalue, $test) = split("\t", $line);
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
