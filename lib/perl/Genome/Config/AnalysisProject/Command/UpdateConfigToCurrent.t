#! /gsc/bin/perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Config::AnalysisProject::Command::UpdateConfigToCurrent') or die;

my $ap = Genome::Config::AnalysisProject->__define__(
    name => '__TEST_AnP__',
    status => 'In Progress',
);
ok($ap, 'define AnP');
my $menu_item = Genome::Config::AnalysisMenu::Item->__define__(
    name => '__DEFAULT__',
);
ok($menu_item, 'define menu item');
my $profile_item = Genome::Config::Profile::Item->__define__(
    status => 'active',
    analysis_project => $ap,
    analysis_menu_item => $menu_item,
);
ok($profile_item, 'define config item');

# Non contrete profile item - no op
my $cmd = Genome::Config::AnalysisProject::Command::UpdateConfigToCurrent->create(profile_item => $profile_item);
ok($cmd, 'create command');
ok($cmd->execute, 'execute command');
is($cmd->status_message, 'This profile item has not yet been concretized, replacing it with the current defaults is a no-op.', 'correct status msg');

my @profile_items = sort { $a->created_at cmp $b->created_at } $ap->config_items;
is(@profile_items, 1, 'AnP has one profile item');
is($profile_items[0]->id, $profile_item->id, 'AnP has the old profile item');
is($profile_items[0]->status, 'active', 'old profile item status is active');

# Contrete profile item - creates new one
my $allocation = Genome::Disk::Allocation->__define__(owner => $profile_item);
$cmd = Genome::Config::AnalysisProject::Command::UpdateConfigToCurrent->create(profile_item => $profile_item);
ok($cmd, 'create command');
ok($cmd->execute, 'execute command');
like($cmd->status_message, qr/Successfully updated/, 'correct status msg');

@profile_items = sort { $a->created_at cmp $b->created_at } $ap->config_items;
is(@profile_items, 2, 'AnP has two profile items');
is($profile_items[0]->id, $profile_item->id, 'AnP has old profile item');
isnt($profile_items[1]->id, $profile_item->id, 'Anp has new profile item');
is($profile_items[0]->status, 'disabled', 'old profile item status is disabled');
is($profile_items[1]->status, 'active', 'new profile item status is active');
is($profile_items[0]->analysis_menu_item, $profile_items[1]->analysis_menu_item, 'profile items have the same menu item');

done_testing();
