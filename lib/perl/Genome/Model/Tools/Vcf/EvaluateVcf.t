#!/usr/bin/env genome-perl

# P R A G M A S ###############################################################
use strict;
use warnings;

# M O D U L E S ###############################################################
use Test::More;
use Path::Class;
use above 'Genome';
#use Genome::Model::Tools::Vcf::EvaluateVcf;

# M A I N #####################################################################

#This is a bare minimum test that just compiles Perl and the UR class.
use_ok('Genome::Model::Tools::Vcf::EvaluateVcf');

SKIP: {
    skip "this is a long-running test", 18;
    run_evaluate_vcf();
}

done_testing();

# S U B R O U T I N E S #######################################################
sub run_evaluate_vcf {
    my ($vcf, $gold_vcf, $roi, $true_negative_bed, $output_dir, $ref) =
      setup_evaluate_vcf_params();

    my $cmd = Genome::Model::Tools::Vcf::EvaluateVcf->create(
        vcf               => $vcf->stringify,
        old_sample        => "H_IJ-NA12878",
        new_sample        => "NA12878",
        gold_vcf          => $gold_vcf->stringify,
        gold_sample       => "NA12878",
        roi               => $roi->stringify,
        true_negative_bed => $true_negative_bed->stringify,
        output_directory  => $output_dir->stringify,
        reference         => $ref,
    );
    ok($cmd, 'Got a Genome::Model::Tools::Vcf::EvaluateVcf instance');
    diag('executing command');
    ok($cmd->execute, 'successfully ran execute');
    my $stats = $cmd->rawstats;
    ok($stats, 'got the raw statistics from vcf evaluation');

    is($cmd->stat_true_positive_found_exact(),
        28820, 'true_positive_found_exact passes');

    is($cmd->stat_total_true_positive_exact(),
        2741358, 'total_true_positive_exact passes');

    is($cmd->stat_sensitivity_exact(),
        0.0105130376988339, 'sensitivity_exact passes');

    is($cmd->stat_true_positive_found_partial(),
        28843, 'true_positive_found_partial passes');

    is($cmd->stat_total_true_positive_partial(),
        2741984, 'total_true_positive_partial passes');

    is($cmd->stat_sensitivity_partial(),
        0.0105190256398287, 'sensitivity_partial passes');

    is($cmd->stat_false_positive_exact(), 139, 'false_positive_exact passes');

    is($cmd->stat_false_positive_partial(),
        116, 'false_positive_partial passes');

    is($cmd->stat_true_negatives(), 2191931600, 'true_negatives passes');

    is($cmd->stat_exact_specificity(),
        0.999999936585612, 'exact_specificity passes');

    is($cmd->stat_partial_specificity(),
        0.999999947078641, 'partial_specificity passes');

    is($cmd->stat_exact_ppv(),   0.995200110501053, 'exact_ppv passes');

    is($cmd->stat_partial_ppv(), 0.995994336821023, 'partial_ppv passes');
    
    is($cmd->stat_vcf_lines_overlapping_tn(),
        47, 'vcf_lines_overlapping_tn passes');

    is($cmd->stat_lines_specificity_in_tn_only(),
        0.999999978557725, 'lines_specificity_in_tn_only passes');
}

sub setup_evaluate_vcf_params {
    # inputs
    my $base_dir =
      Path::Class::Dir->new('/gscmnt/gc2801/analytics/tabbott/vcf-evaluate');

    # Input VCF is based on original vcf file:
    #/gscmnt/gc13028/info/model_data/6195f53ff81046959a3bfee124c186f6/buildf60b23fdb26c45f486cf2d5937a7502b/variants/snvs.vcf.gz";
    my $vcf = $base_dir->file("src/vcf-evaluation/tmp/small.vcf.gz");

    my $gold_vcf = $base_dir->file("gold/snvs.vcf.gz");
    my $roi = $base_dir->file("gold/highconf.bed.gz");
    my $true_negative_bed = $base_dir->file("gold/tn.bed.gz");

    # output directory
    my $output_dir = Path::Class::Dir->new('/tmp/vcf/test');
    if (-e $output_dir) {
        $output_dir->rmtree;
    }
    $output_dir->mkpath;

    my $ref = Genome::Model::Build::ReferenceSequence->get(106942997);

    return ($vcf, $gold_vcf, $roi, $true_negative_bed, $output_dir, $ref);
}
