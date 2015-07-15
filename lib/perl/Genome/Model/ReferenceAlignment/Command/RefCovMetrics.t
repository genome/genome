#!/usr/bin/env genome-perl

use strict;
use warnings 'FATAL';

use above 'Genome';

require Test::MockObject;
use Test::More;
plan tests => 4;

my $class = 'Genome::Model::ReferenceAlignment::Command::RefCovMetrics';
use_ok($class) or die;

my $normal_result = _create_refcov_result_for_normal();
my $tumor_result = _create_refcov_result_for_tumor();

my $build = Test::MockObject->new();
$build->set_always('is_current', 1);

my $model = Test::MockObject->new();
$model->set_isa('Genome::Model::ReferenceAlignment', 'Genome::Model');
$model->set_always('name', 'JohnnyVegas');
$model->set_always('__display_name__', 'JohnnyVegas (1)');
$model->set_always('subject', Genome::Sample->__define__(name => '_TEST_SAMPLE_', common_name => 'tumor'));
$model->set_list('builds', $build);

subtest 'no results' => sub{
    plan tests => 2;

    my $cmd = $class->execute(models => $model, type => 'coverage');
    ok($cmd->result, 'execute for coverage');
    is($cmd->warning_message, 'No results for model: '.$model->__display_name__, 'correct warning mesage');

};

subtest 'coverage from ref-align' => sub{
    plan tests => 2;

    $build->set_list('results', $normal_result);

    my $out = Genome::Sys->create_temp_file_path('model.coverage');
    my $cmd = $class->execute(models => $model, type => 'coverage', output_path => $out);
    ok($cmd->result, 'execute for coverage');
    my @lines = Genome::Sys->read_file($out);
    is_deeply(
        \@lines,
        [
            join("\t", (qw/ model_name coverage-wingspan_0_1_mean_depth coverage-wingspan_0_1_target_base_pair coverage-wingspan_0_40_mean_depth coverage-wingspan_0_40_target_base_pair coverage-wingspan_500_1_mean_depth coverage-wingspan_500_1_target_base_pair coverage-wingspan_500_40_mean_depth coverage-wingspan_500_40_target_base_pair /))."\n",
            join("\t", (qw/ JohnnyVegas 58.303 368332 46.158 368332 17.683 2302332 13.115 2302332 /))."\n",
        ],
        'output matches',
    );

};

subtest 'alignment from som-val' => sub{
    plan tests => 2;

    $model->set_isa('Genome::Model::SomaticValidation', 'Genome::Model');
    $build->set_list('results', $normal_result, $tumor_result);

    my $out = Genome::Sys->create_temp_file_path('model.alignment');
    my $cmd = $class->execute(models => $model, type => 'alignment', row_ids => [qw/ sample_name sample_common_name result_id /], output_path => $out);
    ok($cmd->result, 'execute for alignment');
    my @lines = Genome::Sys->read_file($out);
    is_deeply(
        \@lines,
        [
        join("\t", (qw/ sample_name	sample_common_name result_id alignment-wingspan_0_total_bp	alignment-wingspan_0_unique_off_target_aligned_bp	alignment-wingspan_0_unique_target_aligned_bp	alignment-wingspan_500_total_bp	alignment-wingspan_500_unique_off_target_aligned_bp	alignment-wingspan_500_unique_target_aligned_bp /))."\n",
        join("\t", (qw/ SAMPLE-N normal 11 11 11 11 11 11 11 /))."\n",
        join("\t", (qw/ SAMPLE-T tumor 22 22 22 22 22 22 22 /))."\n",
        ],
        'output matches',
    );

};

done_testing();

###

sub _alignment_summary_hash_ref {
    my $id = shift;
    return {
        '0' => {
            'total_bp' => $id,
            'unique_off_target_aligned_bp' => $id,
            'unique_target_aligned_bp' => $id,
        },
        '500' => {
            'total_bp' => $id,
            'unique_off_target_aligned_bp' => $id,
            'unique_target_aligned_bp' => $id,
        }
    };
}

sub _create_refcov_result_for_normal {
    my $sample = Test::MockObject->new();
    $sample->set_always('name', 'SAMPLE-N');
    $sample->set_always('common_name', 'normal');

    my $instrument_data = Test::MockObject->new;
    $instrument_data->set_always('sample', $sample);

    my $ar = Test::MockObject->new();
    $ar->set_always('instrument_data', $instrument_data);

    my $result = Test::MockObject->new();
    $result->set_isa('Genome::InstrumentData::AlignmentResult::Merged::CoverageStats');
    $result->set_always('id', 11);
    $result->set_always('alignment_result', $ar);
    $result->set_always('coverage_stats_summary_hash_ref', {
            0 => {
                1  => { target_base_pair => 368332, mean_depth => 58.303, },
                40 => { target_base_pair => 368332, mean_depth => 46.158, },
            },
            500 => {
                1  => { target_base_pair => 2302332, mean_depth => 17.683, },
                40 => { target_base_pair => 2302332, mean_depth => 13.115, },

            }
        },
    );
    $result->set_always('alignment_summary_hash_ref', _alignment_summary_hash_ref($result->id));

    return $result;
}

sub _create_refcov_result_for_tumor {
    my $sample = Test::MockObject->new();
    $sample->set_always('name', 'SAMPLE-T');
    $sample->set_always('common_name', 'tumor');

    my @instrument_data = map { my $o = Test::MockObject->new; $o->set_always('sample', $sample); $o } ( 1, 2 );

    my $ar = Test::MockObject->new();
    $ar->set_always('instrument_data', @instrument_data);

    my $result = Test::MockObject->new();
    $result->set_isa('Genome::InstrumentData::AlignmentResult::Merged::CoverageStats');
    $result->set_always('id', 22);
    $result->set_always('alignment_result', $ar);
    $result->set_always('alignment_summary_hash_ref', _alignment_summary_hash_ref($result->id));

    return $result;
}

