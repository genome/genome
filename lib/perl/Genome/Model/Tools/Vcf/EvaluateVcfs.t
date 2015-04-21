#!/usr/bin/env genome-perl

# P R A G M A S ###############################################################
use strict;
use warnings;

# M O D U L E S ###############################################################
use Test::More;
use Test::Number::Delta within => 1e-4;
use Path::Class;
use YAML::XS;
use above 'Genome';

# M A I N #####################################################################

#This is a bare minimum test that just compiles Perl and the UR class.
use_ok('Genome::Model::Tools::Vcf::EvaluateVcfs');

# This directory contains the test data set (based from BIO-1176)
my $basedir = Path::Class::Dir->new(
    '/gscmnt/gc2801/analytics/idas',
    'jira/BIO-1387/cle-data'
);

SKIP: {
    skip "this is a long-running test & must be run within TGI", 75;

    my $cmd = run_evaluate_vcfs();
    check_stats($cmd);
}

done_testing();

# S U B R O U T I N E S #######################################################
sub run_evaluate_vcfs {
    my %params = setup_evaluation_params();

    my $cmd = Genome::Model::Tools::Vcf::EvaluateVcfs->create(%params);
    ok($cmd, 'Got a Genome::Model::Tools::Vcf::EvaluateVcfs instance');
    diag('executing command');
    ok($cmd->execute, 'successfully ran execute');
    return $cmd;
}

sub check_stats {
    my $cmd = shift;

    my $expected_results = get_expected_results();
    ok($expected_results, "Loaded the expected results");

    my $rawdata = $cmd->rawdata();
    ok($rawdata, 'Got the raw output command data');

    my @stat_names = $cmd->stat_types;

    for my $set (@{$rawdata}) {
        my ($name, $type) = ($set->{'name'}, $set->{'variant_type'});
        diag("\nvalidating stat set -- name: '$name' | type: '$type'\n\n");

        my ($expected_set) =
          grep { $_->{'name'} eq $name && $_->{'variant_type'} eq $type }
          @{$expected_results};
        ok($expected_set, "Got expected results set");

        for my $stat (@stat_names) {
            my $calculated = $set->{'stats'}->{$stat};
            my $expected = $expected_set->{'stats'}->{$stat};
            delta_ok($calculated, $expected, "stat: '$stat' matches up");
        }
    }
}

sub setup_evaluation_params {
    my $tmp = Path::Class::Dir->new(Genome::Sys->create_temp_directory);

    # inputs

    my $outdir = setup_output_dir($tmp);
    ok (-e $outdir, "Created an output directory : $outdir");

    my $config = $basedir->file('edited_germline_config.test.tsv');
    ok(-e "$config", "Got the config file : $config");

    my $roi = $basedir->file('NGv3_RMG1_intersect_NIST.bed');
    ok(-e "$roi", "Got the roi file : $roi");

    my $gold_snvs_vcf = $basedir->file('gold-snvs.vcf.gz');
    ok(-e "$gold_snvs_vcf", "Got the gold snvs vcf file : $gold_snvs_vcf");

    my $gold_indels_vcf = $basedir->file('gold-indels.vcf.gz');
    ok(-e "$gold_indels_vcf", "Got the gold indels vcf file : $gold_indels_vcf");

    my $true_negatives_bed = $basedir->file('true-negatives.bed.gz');
    ok(-e "$true_negatives_bed", "Got the true negatives bed file : $true_negatives_bed");

    my %params = (
        output_directory       => "$outdir",
        config_file            => "$config",
        roi                    => "$roi",
        gold_snv_vcf           => "$gold_snvs_vcf",
        gold_indel_vcf         => "$gold_indels_vcf",
        true_negative_bed      => "$true_negatives_bed",
        gold_sample            => 'NA12878',
        pass_filter_expression => "-f 'ALLFILTERSPASS > 0'",
#        true_negative_size     => 45235743,
    );

    return %params;
}

sub get_expected_results {
    my $expected = $basedir->file('answers.yml');
    ok(-e "$expected", "Got the expected results file : $expected");

    my $results = YAML::XS::LoadFile($expected);
    return $results;
}

sub setup_output_dir {
    my $tmp = shift;
    my $outdir = $tmp->subdir('evals');
    $outdir->mkpath;
    return $outdir;
}

__END__
