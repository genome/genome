#!/usr/bin/env genome-perl
use strict;
use warnings;
use above "Genome";
use Test::More;

no warnings 'redefine';
*Genome::Sys::current_user_is_admin = sub {return 1 };
use warnings;

use_ok("Genome::Sys::User") or die "cannot contiue w/o the user module";

my $u0 = Genome::Sys::User->create(
    email => 'someone@somewhere.org',
    name => 'Fake McFakerston'
);
ok($u0,'create a user object');
is($u0->id, 'someone@somewhere.org', 'user id is email');

my $u1 = Genome::Sys::User->get(email => 'someone@somewhere.org');
ok($u1, 'got a user object');
is_deeply($u1, $u0, 'get user');

my $project = Genome::Project->create(name => 'user_test_project');
$u0->add_project($project);

is_deeply([$u0->projects], [$project], "Got project from project creator");

done_testing();
