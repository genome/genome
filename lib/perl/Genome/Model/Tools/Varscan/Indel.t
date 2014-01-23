#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use File::Spec;

use Genome::Model::Tools::TestHelpers::General qw(
    get_test_dir
    compare_to_blessed_file
);
use Genome::Model::Tools::TestHelpers::Data qw(
    get_shared_test_data
);

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

my $class = 'Genome::Model::Tools::Varscan::Indel';
use_ok($class);

my $VERSION = 1;
my $SHARED_DATA_VERSION = 1;

my $ref_fasta = get_ref_fasta();

my $test_dir = get_test_dir($class, $VERSION);
my $input_bam = File::Spec->join($test_dir, '134053361-region-limited.bam');
ok(-s $input_bam, "Found input Bam file: ($input_bam)") || die;

my $output_directory = Genome::Sys->create_temp_directory();

test_varscan_command();
test_without_output_vcf();
test_with_output_vcf();

done_testing();

sub get_ref_fasta {
    return get_shared_test_data('human_reference_37.fa', $SHARED_DATA_VERSION);
}

sub test_without_output_vcf {
    my $output_file = File::Spec->join($output_directory, 'output.cns');
    my $cmd = $class->create(
        bam_file => $input_bam,
        ref_fasta => $ref_fasta,
        output_file => $output_file,
    );
    ok($cmd->execute(), 'Successfully, ran command');
    compare_to_blessed_file(
        output_path => $cmd->output_file,
        output_dir => $output_directory,
        test_dir => $test_dir,
    );
}

sub test_with_output_vcf {
    my $output_file = File::Spec->join($output_directory, 'output.vcf');
    my $cmd = $class->create(
        bam_file => $input_bam,
        ref_fasta => $ref_fasta,
        output_file => $output_file,
        output_vcf => 1,
        vcf_sample_name => 'vcf_sample_name',
    );
    ok($cmd->execute(), 'Successfully, ran command');
    compare_to_blessed_file(
        output_path => $cmd->output_file,
        output_dir => $output_directory,
        test_dir => $test_dir,
    );
}

sub get_command {
    my %additional_params = @_;
    my %params = (
        bam_file => 'bam_file',
        output_file => 'output_file',
        ref_fasta => 'ref_fasta',
    );
    return $class->create(%params, %additional_params);
}

sub test_varscan_command {
    my $cmd = get_command(output_vcf => 0);
    is('mpileup2indel', $cmd->varscan_command, 'Correctly uses mpileup2indel if output_vcf => 0');

    $cmd = get_command(output_vcf => 1);
    is('mpileup2indel', $cmd->varscan_command, 'Correctly uses mpileup2indel if output_vcf => 1');
}

1;
