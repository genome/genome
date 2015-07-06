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
$model->set_always('name', 'Johnny Vegas');
$model->set_always('builds', ($build));

subtest 'coverage' => sub{
    plan tests => 1;
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

    my $cmd = $class->execute(models => $model, type => 'coverage');
    ok($cmd->result, 'execute for coverage');


};

subtest 'alignment' => sub{
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
    my $cmd = $class->execute(models => $model, type => 'alignment');
    ok($cmd->result, 'execute for alignment');

};

done_testing();
