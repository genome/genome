#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Data::Dumper 'Dumper';
use Test::More;

use_ok('Genome::Command::Crud') or die;
use_ok('Genome::Command::Create') or die;
use_ok('Genome::Command::Update') or die;
use_ok('Genome::Command::Delete') or die;

# class

class SuperHero {
    is => 'UR::Object',
    has => [
        name => { is => 'Text', doc => 'Name of the person', },
        special_power => { is => 'Text', doc => 'Superpower', },
        ],
};
sub SuperHero::__display_name__ {
    return $_[0]->name;
}

# INIT
my %config = (
    target_class => 'SuperHero',
    update => { },
);

ok(Genome::Command::Crud->init_sub_commands(%config), 'init crud commands') or die;
# MAIN TREE
my $main_tree_meta = SuperHero::Command->__meta__;
ok($main_tree_meta, 'MAIN TREE meta');

# CREATE
# meta 
my $create_meta = SuperHero::Command::Create->__meta__;
ok($create_meta, 'CREATE meta');

is(SuperHero::Command::Create->_target_name, 'super hero', 'CREATE: _target_name');
is(SuperHero::Command::Create->_target_class, 'SuperHero', 'CREATE: _target_class');

# Create Spiderman
my $create_spiderman = SuperHero::Command::Create->create(
    name => 'Spiderman',
    special_power => 'Shoots webs',
);

ok($create_spiderman, "CREATE command created");
#my $spiderman = $create_spiderman->execute();
my $spiderman = SuperHero->create(
    name => 'Spiderman',
    special_power => 'Shoots webs',
);
ok($spiderman, "Spiderman created");
is($spiderman->name, 'Spiderman', "Spiderman's name is correct");

# LIST - this is in UR::Object::Command...
my $list_meta = SuperHero::Command::List->__meta__;
ok($list_meta, 'LIST meta');

# UPDATE
# meta
my $update_meta = SuperHero::Command::Update->__meta__;
ok($update_meta, 'update meta');

is(SuperHero::Command::Update::Name->_target_name_pl, 'super heroes', 'UPDATE: _target_name_pl');
is(SuperHero::Command::Update::Name->_target_name_pl_ub, 'super_heroes', 'UPDATE: _target_name_pl_ub');

# update name
is($spiderman->name, 'Spiderman', 'Spiderman has a name');
my $update_success = SuperHero::Command::Update::Name->create(
    super_heroes => [ $spiderman ],
    value => 'Spidey',
);
ok($update_success, 'UPDATE create Change Spideys name');
$update_success->dump_status_messages(1);
ok($update_success->execute, "UPDATE execute");
is($spiderman->name, 'Spidey', 'Spiderman has a new name');

# DELETE
# meta
my $delete_meta = SuperHero::Command::Delete->__meta__;
ok($delete_meta, 'DELETE meta');

is(SuperHero::Command::Delete->_target_name_pl, 'super heroes', 'DELETE: _target_name_pl');
is(SuperHero::Command::Delete->_target_name_pl_ub, 'super_heroes', 'DELETE: _target_name_pl_ub');

# success
my $delete_success = SuperHero::Command::Delete->create(
    super_heroes => [ $spiderman ],
);
ok($delete_success, 'DELETE create: Spidey');
$delete_success->dump_status_messages(1);
ok($delete_success->execute, "DELETE execute");
my $deleted_jimmy = SuperHero->get(name => 'Spidey');
ok(!$deleted_jimmy, 'deleted Spiderman confirmed');

# COMMIT
ok(UR::Context->commit, 'commit');

done_testing();
exit;

