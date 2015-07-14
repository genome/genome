#!/usr/bin/env genome-perl

use strict;
use warnings 'FATAL';

use above 'Genome';

require Test::MockObject;
use Test::More;
plan tests => 4;

my $class = 'Genome::Model::ReferenceAlignment::Command::RefCovMetrics';
use_ok($class) or die;

my $build = Test::MockObject->new();
$build->set_always('is_current', 1);
$build->set_series('results');

my $model = Test::MockObject->new();
$model->set_isa('Genome::Model::ReferenceAlignment', 'Genome::Model');
$model->set_always('name', 'JohnnyVegas');
$model->set_always('__display_name__', 'JohnnyVegas (1)');
$model->set_always('subject', Genome::Sample->__define__(name => '_TEST_SAMPLE_', common_name => 'tumor'));
$model->set_series('builds', $build);

subtest 'no results' => sub{
    plan tests => 2;

    my $cmd = $class->execute(models => $model, type => 'coverage');
    ok($cmd->result, 'execute for coverage');
    is($cmd->warning_message, 'No results for model: '.$model->__display_name__, 'correct warning mesage');

};

my $result = Test::MockObject->new();
$result->set_isa('Genome::InstrumentData::AlignmentResult::Merged::CoverageStats');
$result->set_always('id', 11);
$build->set_series('results', $result);

my $instrument_data = Test::MockObject->new;
$instrument_data->set_always('sample_name', 'SAMPLE'.$_);
$instrument_data->set_always('common_name', 'normal');
$result->set_series('instrument_data', $instrument_data);

subtest 'coverage from ref-align' => sub{
    plan tests => 2;

    my $coverage_stats_summary = {
        0 => {
            1  => { target_base_pair => 368332, mean_depth => 58.303, },
            40 => { target_base_pair => 368332, mean_depth => 46.158, },
        },
        500 => {
            1  => { target_base_pair => 2302332, mean_depth => 17.683, },
            40 => { target_base_pair => 2302332, mean_depth => 13.115, },

        }
    };
    $result->set_always('coverage_stats_summary_hash_ref', $coverage_stats_summary);

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

$model->set_isa('Genome::Model::SomaticValidation', 'Genome::Model');
my $tumor_result = Test::MockObject->new();
$tumor_result->set_isa('Genome::InstrumentData::AlignmentResult::Merged::CoverageStats');
$tumor_result->set_always('id', 22);
$build->set_series('results', $result, $tumor_result);

my @instrument_data;
for (0..1) {
    push @instrument_data, Test::MockObject->new;
    $instrument_data[$_]->set_always('sample_name', 'SAMPLE-T');
    $instrument_data[$_]->set_always('common_name', 'tumor');
}
$result->set_series('instrument_data', @instrument_data);

subtest 'alignment from som-val' => sub{
    plan tests => 2;

    my $alignment_summary_hash_ref = {
        '0' => {
            'total_bp' => '123456789',
            'unique_off_target_aligned_bp' => '12345',
            'unique_target_aligned_bp' => '6789',
        },
        '500' => {
            'total_bp' => '987654321',
            'unique_off_target_aligned_bp' => '98765',
            'unique_target_aligned_bp' => '4321',
        }
    };
    $tumor_result->set_always('alignment_summary_hash_ref', $alignment_summary_hash_ref);

    my $out = Genome::Sys->create_temp_file_path('model.alignment');
    my $cmd = $class->execute(models => $model, type => 'alignment', row_ids => [qw/ sample_name sample_common_name result_id /], output_path => $out);
    ok($cmd->result, 'execute for alignment');
    my @lines = Genome::Sys->read_file($out);
    is_deeply(
        \@lines,
        [
            join("\t", (qw/ sample_name	sample_common_name result_id alignment-wingspan_0_total_bp	alignment-wingspan_0_unique_off_target_aligned_bp	alignment-wingspan_0_unique_target_aligned_bp	alignment-wingspan_500_total_bp	alignment-wingspan_500_unique_off_target_aligned_bp	alignment-wingspan_500_unique_target_aligned_bp /))."\n",
            join("\t", (qw/ SAMPLE-N normal 11 366769868 299651882 29092686 366769868 298629738 30114830 /))."\n",
            join("\t", (qw/ SAMPLE-T tumor 22 123456789 12345 6789 987654321 98765 4321 /))."\n",
        ],
        'output matches',
    );

};

done_testing();
