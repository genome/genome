#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Genome::Sys;

use Test::More;

use_ok('Genome::Model::Tools::Gatk::WithNumberOfThreads') or die;

class GMT::Gatk::TestWithNumberThreads {
    is => 'Genome::Model::Tools::Gatk::WithNumberOfThreads',
};
sub GMT::Gatk::TestWithNumberThreads::_build_gatk_command { return 'gatk' };

my $cmd = GMT::Gatk::TestWithNumberThreads->create();
is($cmd->build_gatk_command, 'gatk', 'correct build_gatk_command w/o threads');
$cmd->number_of_threads(8);
is($cmd->build_gatk_command, 'gatk -nt 8', 'correct build_gatk_command w/ 8 threads');

# Test shellcmd OK
Sub::Install::reinstall_sub({
        code => sub{ return 1; },
        into => "Genome::Sys",
        as   => 'shellcmd',
    });
my $rv = $cmd->execute;
ok($rv, '_add_threads_option_to_command_and_run succeeds');
is($cmd->shellcmd_exit_code, 0, 'correct shellcmd_exit_code after successful execute');

# Test shellcmd fail
$cmd = GMT::Gatk::TestWithNumberThreads->create();
Sub::Install::reinstall_sub({
        code => sub{ Carp::croak("ERROR RUNNING COMMAND.  Exit code 134 from: cmd\nSee the command's captured STDERR (if it exists) for more information"); },
        into => "Genome::Sys",
        as   => 'shellcmd',
    });
$rv = $cmd->execute;
ok(!$rv, 'execute fails');
is($cmd->shellcmd_exit_code, 134, 'correct shellcmd_exit_code after failed execute');

done_testing();
