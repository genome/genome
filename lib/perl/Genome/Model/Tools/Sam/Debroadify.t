#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
use File::Compare;
use File::Temp;
use above 'Genome';

my $class = 'Genome::Model::Tools::Sam::Debroadify';
use_ok($class);

my %test_env = setup_test_env();
test_sam_inputs();
test_reference_input();
test_debroadify_sam();
test_debroadify_bam();

done_testing();

sub test_sam_inputs {
    do {
        my $cmd = $class->create(
            input_file => $test_env{input_sam_file} . '',
            output_file => $test_env{output_dir} . '/output.sam',
        );
        my $inputs_did_validate = eval { $cmd->validate_inputs() };
        is($inputs_did_validate, 1, 'valid input sam passes');
    };

    do {
        my $cmd = $class->create(
            input_file => $test_env{output_dir} . '/nonexistant.sam',
            output_file => $test_env{output_dir} . '/output.sam',
        );
        my $inputs_did_validate = eval { $cmd->validate_inputs() };
        is($inputs_did_validate, undef, 'invalid input sam fails');
    };
}

sub test_reference_input {
    do {
        my $cmd = $class->create(
            input_file => $test_env{input_bam_file} . '',
            output_file => $test_env{output_dir} . '/output.bam',
            reference_file => $test_env{reference_file} . '',
        );
        my $inputs_did_validate = eval { $cmd->validate_inputs() };
        is($inputs_did_validate, 1, 'valid reference passes');
    };

    do {
        my $cmd = $class->create(
            input_file => $test_env{input_bam_file} . '',
            output_file => $test_env{output_dir} . '/output.bam',
            reference_file => $test_env{output_dir} . '/nonexistant.fa',
        );
        my $inputs_did_validate = eval { $cmd->validate_inputs() };
        is($inputs_did_validate, undef, 'invalid reference passes');
    };

    do {
        my $cmd = $class->create(
            input_file => $test_env{input_bam_file} . '',
            output_file => $test_env{output_dir} . '/output.bam',
        );
        my $inputs_did_validate = eval { $cmd->validate_inputs() };
        is($inputs_did_validate, undef, 'missing reference fails with bam output');
    };
}

sub test_debroadify_sam {
    ok(-d $test_env{output_dir}, 'output_dir (' . $test_env{output_dir} . ') exists');

    my $input_bam = $test_env{input_dir} . '/alignment.bam';
    ok(-e $input_bam, 'input_bam exists') || return;

    my $expected_output_sam = $test_env{input_dir} . '/alignment.sam';
    ok(-e $expected_output_sam, 'expected_output_sam exists') || return;

    my $output_sam = $test_env{output_dir} . '/alignment.sam';
    ok(!-e $output_sam, 'output_sam does not exist') || return;

    my $debroadify_cmd = Genome::Model::Tools::Sam::Debroadify->create(
        input_file => $input_bam,
        output_file => $output_sam,
    );
    ok($debroadify_cmd->execute(), 'debroadify command completed successfully') || return;

    my $compare_cmd = Genome::Model::Tools::Sam::Compare->create(
        file1 => $output_sam,
        file2 => $expected_output_sam,
    );
    ok($compare_cmd->execute(), 'output sam matched expected output') || return;

    my $compare_orig_cmd = Genome::Model::Tools::Sam::Compare->create(
        file1 => $output_sam,
        file2 => $input_bam,
    );
    ok(!$compare_orig_cmd->execute(), 'output sam does not match input') || return;
}

sub test_debroadify_bam {
    ok(-d $test_env{output_dir}, 'output_dir (' . $test_env{output_dir} . ') exists');

    my $input_bam = $test_env{input_dir} . '/alignment.bam';
    ok(-e $input_bam, 'input_bam exists') || return;

    my $expected_output_bam = $test_env{input_dir} . '/alignment_rev.bam';
    ok(-e $expected_output_bam, 'expected_output_bam exists') || return;

    my $output_bam = $test_env{output_dir} . '/alignment_rev.bam';
    ok(!-e $output_bam, 'output_bam does not exist') || return;

    my $debroadify_cmd = Genome::Model::Tools::Sam::Debroadify->create(
        input_file => $input_bam,
        output_file => $output_bam,
        reference_file => $test_env{input_dir} . '/all_sequences_rev.fa',
    );
    ok($debroadify_cmd->execute(), 'debroadify command completed successfully') || return;

    my $compare_cmd = Genome::Model::Tools::Sam::Compare->create(
        file1 => $output_bam,
        file2 => $expected_output_bam,
    );
    ok($compare_cmd->execute(), 'output bam matched expected output') || return;

    my $compare_orig_cmd = Genome::Model::Tools::Sam::Compare->create(
        file1 => $output_bam,
        file2 => $input_bam,
    );
    ok(!$compare_orig_cmd->execute(), 'output bam does not match input') || return;
}

sub setup_test_env {
    my %test_env;

    $test_env{input_dir} = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sam-Debroadify';
    $test_env{output_dir} = File::Temp->newdir();

    $test_env{input_sam_file} = File::Temp->new(
        SUFFIX => '.sam',
    );
    $test_env{input_bam_file} = File::Temp->new(
        SUFFIX => '.bam',
    );
    $test_env{reference_file} = File::Temp->new(
        SUFFIX => '.fa',
    );

    return %test_env;
}
