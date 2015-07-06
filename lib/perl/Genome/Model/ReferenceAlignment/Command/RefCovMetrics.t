#!/usr/bin/env genome-perl

use strict;
use warnings 'FATAL';

use above 'Genome';

require Test::MockObject;
use Test::More;

my $class = 'Genome::Model::ReferenceAlignment::Command::RefCovMetrics';
use_ok($class) or die;

my $build = Test::MockObject->new();
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

my $model = Test::MockObject->new();
$model->set_isa('Genome::Model::ReferenceAlignment', 'Genome::Model');
$model->set_always('builds', ($build));

my $cmd = $class->execute(models => $model);
ok($cmd->result, 'execute');

done_testing();
