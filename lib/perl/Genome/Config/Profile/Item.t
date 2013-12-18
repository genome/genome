#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

my $class = 'Genome::Config::Profile::Item';
use_ok($class);

my $menu_item_contents = 'menu item contents';
my $menu_item_file = _create_file_with_contents($menu_item_contents);
my $menu_item = Genome::Config::AnalysisMenu::Item->create(
    file_path => $menu_item_file,
    name => 'test menu item'
);

my $analysis_project = Genome::Config::AnalysisProject->create(name => 'Test Project');

my $profile_item_from_menu_item = Genome::Config::Profile::Item->create(
    analysis_project => $analysis_project,
    analysis_menu_item => $menu_item,
);

ok(!$profile_item_from_menu_item->_is_concrete, 'it shouldnt start as concrete when made from a menu item');
is($menu_item_file, $profile_item_from_menu_item->file_path, 'it should point to the original menu item file');
is($profile_item_from_menu_item->allocation, undef, 'it should not have an allocation yet');

$profile_item_from_menu_item->concretize();

ok($profile_item_from_menu_item->_is_concrete, 'it should now be concretized');
ok($profile_item_from_menu_item->allocation, 'it should now have an allocation');
ok($menu_item_file ne $profile_item_from_menu_item->file_path, 'it should no longer point to the menu item');
is(Genome::Sys->read_file($profile_item_from_menu_item->file_path), $menu_item_contents, 'the new file should match the menu item');


my $custom_config_file_contents = 'custom config - not from menu item';
my $custom_config_file_path = _create_file_with_contents($custom_config_file_contents);

my $profile_item_from_file = Genome::Config::Profile::Item->create_from_file_path(
    analysis_project => $analysis_project,
    file_path => $custom_config_file_path,
);

ok($profile_item_from_file->_is_concrete, 'it should start out as concrete');
ok($profile_item_from_file->allocation, 'it should start out with an allocation');
ok($profile_item_from_file->file_path ne $custom_config_file_path, 'it should not point to the original file');
is(Genome::Sys->read_file($profile_item_from_file->file_path), $custom_config_file_contents, 'it should copy the original file to the new location');
ok(!$profile_item_from_file->analysis_menu_item, 'it should not have a menu item selected');


sub _create_file_with_contents {
    my $contents = shift;
    my $file_path = Genome::Sys->create_temp_file_path();

    Genome::Sys->write_file($file_path, $contents);

    return $file_path;
}

done_testing();
