use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use Test::More;
use above "Genome";

my $p = Genome::Project->create(name => 'Genome::Project::Unwatch Test');
my $user = Genome::Sys::User->create(id => 'johndoe@wustl.edu');
$p->add_watcher($user);

my $remove_cmd = Genome::Project::Command::Unwatch->create(project => $p, user => $user);
ok($remove_cmd->execute(), 'executed remove_cmd');
like($remove_cmd->status_message, qr/Removing/, 'got remove message');

my $already_cmd = Genome::Project::Command::Unwatch->create(project => $p, user => $user);
$DB::single = 1;
ok($already_cmd->execute(), 'executed already_cmd');
like($already_cmd->status_message, qr/not watching/, 'got already watching message');

done_testing();
