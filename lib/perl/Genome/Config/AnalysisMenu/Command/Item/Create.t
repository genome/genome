#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use above 'Genome';
use Test::More tests => 6;
use Test::Exception;

my $class = 'Genome::Config::AnalysisMenu::Command::Item::Create';
use_ok($class);

my $name = "TestMenuItemWithCheese";
my $description = "MEAT AND CHEESE FOR THE MEAT AND CHEESE GODS!";
my $tmp_dir = Genome::Sys->create_temp_directory;
my $file_path = File::Spec->join($tmp_dir, 'test_menu_item.yml');
Genome::Sys->write_file($file_path, "2 all beef patties");

local $ENV{'GENOME_ANALYSIS_PROJECT_DEFAULTS'} = '/dev/null';
my $cmd = $class->create(description => $description, file_path => $file_path, name => $name);
ok($cmd, "Successfully created the command");
dies_ok(sub{$cmd->execute}, 'Command fails when file is not in menu directory');

local $ENV{'GENOME_ANALYSIS_PROJECT_DEFAULTS'} = $tmp_dir;
my $cmd2 = $class->create(description => $description, file_path => $file_path, name => $name);
ok($cmd2, "Successfully created the command");
ok($cmd2->execute, 'Successfully executed the command');
ok(Genome::Config::AnalysisMenu::Item->get(file_path => $file_path), 'Successfully found the Menu Item');

1;
