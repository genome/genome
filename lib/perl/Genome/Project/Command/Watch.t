use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use Test::More;
use above "Genome";

use Genome::Sys;
no warnings 'redefine';
*Genome::Sys::current_user_is_admin = sub { return 1 };
use warnings;

my $p = Genome::Project->create(name => 'Genome::Project::Watch Test');
my $user = Genome::Sys::User->create(id => 'johndoe@wustl.edu');

my $add_cmd = Genome::Project::Command::Watch->create(project => $p, user => $user);
ok($add_cmd->execute(), 'executed add_cmd');
like($add_cmd->status_message, qr/Adding/, 'got add message');

my $already_cmd = Genome::Project::Command::Watch->create(project => $p, user => $user);
$DB::single = 1;
ok($already_cmd->execute(), 'executed already_cmd');
like($already_cmd->status_message, qr/already watching/, 'got already watching message');

done_testing();
