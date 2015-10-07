#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Genome::Utility::Test qw(compare_ok abort);
use Test::More tests => 12;

eval {
    # create a class instance
    my $class = 'Genome::Model::Tools::Analysis::Concordance';
    use_ok($class);

    # check test data files
    my $data_dir = Genome::Utility::Test->data_dir($class, 'v3');
    ok(-d $data_dir, "data_dir exists: $data_dir") or abort;

    # check inputs
    my $bam_file_1 = "$data_dir/168_norm_chr21.10k.bam";
    ok(-s $bam_file_1, 'bam file 1 exists: $bam_file_1') or abort;
    my $bam_file_2 = "$data_dir/168_pre_chr21.10k.bam";
    ok(-s $bam_file_2, 'bam file 2 exists: $bam_file_2') or abort;

    my $reference_fasta = "$data_dir/21.fa";
    ok(-s $reference_fasta, 'reference genome exists: $reference_fasta') or abort;

    my $snp_file = "$data_dir/input.snps";
    ok(-s $snp_file, 'SNP file exists: $snp_file') or abort;

    my $result_file = "$data_dir/expected_outfile";
    ok(-s $result_file, 'Result file exists: $result_file') or abort;

    # create and execute command
    my $tmpdir = Genome::Sys->create_temp_directory();
    my $cmd = $class->create(
        bam_file_1   => $bam_file_1,
        bam_file_2   => $bam_file_2,
        snp_file    => $snp_file,
        reference_fasta    => $reference_fasta,
        output_file => "$tmpdir/out_file",
        output_genotypes => "$tmpdir/out_geno",
        bam_readcount_version => 0.6,
    );

    ok($cmd, 'created command') or abort;
    ok($cmd->execute, 'executed command') or abort;

    # check outputs
    ok(-s "$tmpdir/out_file", 'output_file has size');
    my %compare_args = (
        replace => [
            [ qr(\Q$data_dir\E) => 'TEST_INPUTS_DIR' ],
        ],
    );
    compare_ok("$tmpdir/out_file", "$data_dir/expected_outfile", 'output_file matched expected', %compare_args);
    compare_ok("$tmpdir/out_geno", "$data_dir/genotype.output", 'output_file matched expected', %compare_args);
};
