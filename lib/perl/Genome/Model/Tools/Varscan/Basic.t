#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use File::Spec;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

my $class = 'Genome::Model::Tools::Varscan::Basic';
use_ok($class);

test_samtools_command();
test_output_vcf_string();
test_input_files();

done_testing();

sub get_command {
    my %additional_params = @_;
    my %params = (
        bam_file => 'bam_file',
        output_file => 'output_file',
        ref_fasta => 'ref_fasta',
    );
    return $class->create(%params, %additional_params);
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

1;
