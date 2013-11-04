#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use File::Spec;
# probably a better place for these than in Dindel...
use Genome::Model::Tools::Dindel::TestHelpers qw(
    get_test_dir
    get_ref_fasta
    compare_output_to_test_data
);

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

my $class = 'Genome::Model::Tools::Varscan::Consensus';
use_ok($class);

my $VERSION = 2; # Bump this each time test data changes

my $ref_fasta = get_ref_fasta();

my $test_dir = get_test_dir($class, $VERSION);
my $input_bam = File::Spec->join($test_dir, '134053361-region-limited.bam');
ok(-s $input_bam, "Found input Bam file");

my $output_directory = Genome::Sys->create_temp_directory();

test_without_output_vcf();
test_with_output_vcf();

test_varscan_command();
test_samtools_command();
test_output_vcf_string();
test_input_files();

done_testing();

sub test_without_output_vcf {
    my $output_file = File::Spec->join($output_directory, 'output.cns');
    my $cmd = $class->create(
        bam_file => $input_bam,
        ref_fasta => $ref_fasta,
        output_file => $output_file,
    );
    ok($cmd->execute(), 'Successfully, ran command');
    compare_output_to_test_data($cmd->output_file, $output_directory, $test_dir);
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
    compare_output_to_test_data($cmd->output_file, $output_directory, $test_dir);
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
    is('pileup2cns', $cmd->varscan_command, 'Correctly uses pileup2cns if output_vcf => 0');

    $cmd = get_command(output_vcf => 1);
    is('mpileup2cns', $cmd->varscan_command, 'Correctly uses mpileup2cns if output_vcf => 1');
}

sub test_samtools_command {
    my $cmd = get_command();
    my $expected = '/gsc/pkg/bio/samtools/samtools-0.1.16/samtools mpileup -B -f ref_fasta -q 10  bam_file';
    is($expected, $cmd->samtools_command, 'Correct output from samtools_command');

    $cmd = get_command(position_list_file=>'plf');
    $expected = '/gsc/pkg/bio/samtools/samtools-0.1.16/samtools mpileup -B -f ref_fasta -q 10 -l plf bam_file';
    is($expected, $cmd->samtools_command, 'Correct output from samtools_command with position_list_file specified.');
}

sub test_output_vcf_string {
    my $cmd = get_command(output_vcf => 0);
    is('', $cmd->output_vcf_string, 'Correct output from output_vcf_string with output_vcf => 0');

    $cmd = get_command(output_vcf => 1);
    is('--output-vcf 1', $cmd->output_vcf_string, 'Correct output from output_vcf_string with output_vcf => 1');

    $cmd = get_command(output_vcf => 1, vcf_sample_name => 'vcf_sample_name');
    is('--output-vcf 1 --vcf-sample-list <(echo "vcf_sample_name")', $cmd->output_vcf_string,
        'Correct output with vcf_sample_name');
}

sub test_input_files {
    my $cmd = get_command();
    ok(scalar(grep {$_ eq 'bam_file'} $cmd->input_files), 'bam_file in input_files');
    ok(scalar(grep {$_ eq 'ref_fasta'} $cmd->input_files), 'ref_fasta in input_files');
    ok(!scalar(grep {$_ eq 'plf'} $cmd->input_files), 'position_list_file not in input_files');

    $cmd = get_command(position_list_file=>'plf');
    ok(scalar(grep {$_ eq 'plf'} $cmd->input_files), 'position_list_file is in input_files if specified');
}


