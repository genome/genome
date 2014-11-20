#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::Exception;
use Test::More;

use_ok('Genome::Model::Tools::CgHub::Base') or die;

class Genome::Model::Tools::CgHub::Success {
    is => 'Genome::Model::Tools::CgHub::Base',
};
sub Genome::Model::Tools::CgHub::Success::_build_command { return 'cgquery --help'; };
sub Genome::Model::Tools::CgHub::Success::_verify_success { return 1; };

my $cmd = Genome::Model::Tools::CgHub::Success->create;
ok($cmd, 'create cg hub test cmd');
ok($cmd->execute, 'execute cmd');

class Genome::Model::Tools::CgHub::FailsInBuildCmd {
    is => 'Genome::Model::Tools::CgHub::Base',
};
sub Genome::Model::Tools::CgHub::FailsInBuildCmd::_build_command { return; };
my $fail = Genome::Model::Tools::CgHub::FailsInBuildCmd->create;
ok($fail, 'create cg hub test cmd');
throws_ok(sub{ $fail->execute; }, qr/Failed to build CG Hub command/, 'execute cmd fails in _build_command');

class Genome::Model::Tools::CgHub::FailsInRunCmd {
    is => 'Genome::Model::Tools::CgHub::Base',
};
sub Genome::Model::Tools::CgHub::FailsInRunCmd::_build_command { return 'echo'; };
sub Genome::Model::Tools::CgHub::FailsInRunCmd::_run_command { return; };
$fail = Genome::Model::Tools::CgHub::FailsInRunCmd->create;
ok($fail, 'create cg hub test cmd');
throws_ok(sub{ $fail->execute; }, qr/Failed to run CG Hub command/, 'execute cmd fails in _run_command');

class Genome::Model::Tools::CgHub::FailsInVerifySuccess {
    is => 'Genome::Model::Tools::CgHub::Base',
};
sub Genome::Model::Tools::CgHub::FailsInVerifySuccess::_build_command { return 'echo'; };
sub Genome::Model::Tools::CgHub::FailsInVerifySuccess::_run_command { return 1; };
sub Genome::Model::Tools::CgHub::FailsInVerifySuccess::_verify_success { return; };
$fail = Genome::Model::Tools::CgHub::FailsInVerifySuccess->create;
ok($fail, 'create cg hub test cmd');
throws_ok(sub{ $fail->execute; }, qr/Ran CG Hub command, but it was determined that it was not successful/, 'execute cmd fails in _verify_success');

done_testing();
