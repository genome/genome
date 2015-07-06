#!/usr/bin/env genome-perl

use strict;
use warnings 'FATAL';

use above 'Genome';

require Test::MockObject;
use Test::More;
plan tests => 3;

my $class = 'Genome::Model::ReferenceAlignment::Command::RefCovMetrics';
use_ok($class) or die;

my $build = Test::MockObject->new();
my $model = Test::MockObject->new();
$model->set_isa('Genome::Model::ReferenceAlignment', 'Genome::Model');
$model->set_always('name', 'JohnnyVegas');
$model->set_always('builds', ($build));

subtest 'coverage' => sub{
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
    $build->set_always('coverage_stats_summary_hash_ref', $coverage_stats_summary);

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

subtest 'alignment' => sub{
    plan tests => 2;

    my $alignment_summary_hash_ref = {
        '0' => {
            'total_bp' => '366769868',
            'unique_off_target_aligned_bp' => '299651882',
            'unique_target_aligned_bp' => '29092686',
        },
        '500' => {
            'total_bp' => '366769868',
            'unique_off_target_aligned_bp' => '298629738',
            'unique_target_aligned_bp' => '30114830',
        }
    };
    $build->set_always('alignment_summary_hash_ref', $alignment_summary_hash_ref);

    my $out = Genome::Sys->create_temp_file_path('model.alignment');
    my $cmd = $class->execute(models => $model, type => 'alignment', output_path => $out);
    ok($cmd->result, 'execute for alignment');
    my @lines = Genome::Sys->read_file($out);
    is_deeply(
        \@lines,
        [
            join("\t", (qw/ model_name	alignment-wingspan_0_total_bp	alignment-wingspan_0_unique_off_target_aligned_bp	alignment-wingspan_0_unique_target_aligned_bp	alignment-wingspan_500_total_bp	alignment-wingspan_500_unique_off_target_aligned_bp	alignment-wingspan_500_unique_target_aligned_bp /))."\n",
            join("\t", (qw/ JohnnyVegas 366769868 299651882 29092686 366769868 298629738 30114830 /))."\n",
        ],
        'output matches',
    );

};

done_testing();
