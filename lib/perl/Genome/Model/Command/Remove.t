#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 14;

use_ok('Genome::Model::Command::Remove') or die;

class Genome::Model::Tester { is => 'Genome::ModelDeprecated', };

my $s = Genome::Sample->create(id => -888, name => 'TEST-' . __FILE__ . "-$$");
ok($s, "made a test sample");

my $p = Genome::ProcessingProfile::Tester->create(
    name => "Tester PP for Testing", 
);
ok($p, "made a test processing profile");

my $m = Genome::Model::Tester->create(
    name => 'Testy McTesterson',
    processing_profile_id => $p->id,
    subject => $s,
);
ok($m, "made a test model");

sub ok_run {
    note("running with params @_");
    my $exit_code1 = eval { Genome::Model::Command::Remove->_execute_with_shell_params_and_return_exit_code(@_); };
    ok(!$@, "the command did not crash") or diag($@);
    is($exit_code1, 0, "command believes it succeeded");
}

sub ok_fail {
    note("running with params @_");
    my $exit_code1 = eval { Genome::Model::Command::Remove->_execute_with_shell_params_and_return_exit_code(@_); };
    ok(!$@, "the command did not crash") or diag($@);
    ok($exit_code1 != 0, "command returned failed exit code $exit_code1 as expected");
}

sub ok_crash {
    note("running with params @_");
    my $exit_code1 = eval { Genome::Model::Command::Remove->_execute_with_shell_params_and_return_exit_code(@_); };
    ok($@, "the command did crash with error $@");
    is($exit_code1, undef, "command exit code is undef as expected") or diag($exit_code1); 
}

# attempting to delete nothing will raise an error
ok_fail();

# disable real deletions, and get control over the return value
my $delete_retval = 1;
my $old_delete = \&Genome::Model::delete;
my $new_delete = sub { return $delete_retval };
*Genome::Model::delete = $new_delete;

# check defaults
my $remove_cmd = Genome::Model::Command::Remove->create(models => [$m]);
ok($remove_cmd,'Model remove command created');
ok(!$remove_cmd->force_delete,'Do not force delete model');
$remove_cmd->force_delete(1);
ok($remove_cmd->force_delete,'force delete model');

# the model delete will cause this to fail
$delete_retval = undef;
my $rv;
eval { $rv = $remove_cmd->execute };
ok(!$rv,'remove model did not work because model delete failed');

# but this will work
$delete_retval = 1;
$remove_cmd = Genome::Model::Command::Remove->create(
                                                     models => [$m],
                                                     force_delete => 1,
                                                 );
ok($remove_cmd->execute,'delete model did work');

# restore the original delete logic
*Genome::Model::delete = $old_delete;

# delete the model made above using the command-line pattern match interface
ok_run($m->name,'--force-delete');
isa_ok($m,"UR::DeletedRef", "model object is deleted");

done_testing();
