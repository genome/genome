#!/usr/bin/env genome-perl
use strict;
use warnings;

use Test::More;
use above "Genome";
use Genome::Utility::Test;

my $class = 'Genome::Config::CopyStrategy::TreeCopy';

use_ok($class);
my $source_data_dir = Genome::Utility::Test->data_dir($class, 1);
ok(-d $source_data_dir, "data_dir exists: $source_data_dir") or die;

my $dest_data_dir = Genome::Sys->create_temp_directory();

sub _run_in_eval {
    my ($file, $source_root, $dest_root, $message) = @_;
    eval {
        Genome::Config::CopyStrategy::TreeCopy->copy_config($file, $source_root, $dest_root);
    };
    ok($@, $message);
}

_run_in_eval('this_file_doesnt_exist', $source_data_dir,
    $dest_data_dir, 'dies when given a nonexistent file');

_run_in_eval($source_data_dir . '/level1/test.yml', 'this_source_dir_doesnt_exist',
    $dest_data_dir, 'dies when given a nonexistent source dir');

_run_in_eval($source_data_dir . '/level1/test.yml', $source_data_dir,
    'this_destination_dir_doesnt_exist', 'dies when given a nonexistent destination dir');


my $file1 = "$source_data_dir/level1/test";
Genome::Config::CopyStrategy::TreeCopy->copy_config($file1, $source_data_dir, $dest_data_dir);
ok(-f "$dest_data_dir/level1/test", "copies file over when one level deep");

my $file2 = "$source_data_dir/level1/level2/level3/test";
Genome::Config::CopyStrategy::TreeCopy->copy_config($file2, $source_data_dir, $dest_data_dir);
ok(-f "$dest_data_dir/level1/level2/level3/test", "copies file over when multiple levels deep");


done_testing();
