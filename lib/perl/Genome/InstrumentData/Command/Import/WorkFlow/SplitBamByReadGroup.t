#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

require Genome::Utility::Test;
require File::Temp;
use Test::More;

my $class = 'Genome::InstrumentData::Command::Import::WorkFlow::SplitBamByReadGroup';
use_ok($class) or die;
my $test_dir = Genome::Utility::Test->data_dir_ok($class, 'v1') or die;

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $one_read_group_basename = 'one-read-group.bam';
my $one_read_group_bam_path = $tmp_dir.'/'.$one_read_group_basename;
Genome::Sys->create_symlink($test_dir.'/'.$one_read_group_basename, $one_read_group_bam_path);
ok(-s $one_read_group_bam_path, 'linked one read group bam');
my $cmd = Genome::InstrumentData::Command::Import::WorkFlow::SplitBamByReadGroup->execute(bam_path => $one_read_group_bam_path);
ok($cmd, 'execute');

my $two_read_groups_basename = 'two-read-groups.bam';
my $two_read_groups_bam_path = $tmp_dir.'/'.$two_read_groups_basename;
Genome::Sys->create_symlink($test_dir.'/'.$two_read_groups_basename, $two_read_groups_bam_path);
Genome::Sys->create_symlink($test_dir.'/'.$two_read_groups_basename.'.flagstat', $two_read_groups_bam_path.'.flagstat');
ok(-s $two_read_groups_bam_path, 'linked two read groups bam');
$cmd = Genome::InstrumentData::Command::Import::WorkFlow::SplitBamByReadGroup->execute(bam_path => $two_read_groups_bam_path);
ok($cmd, 'execute');

<STDIN>;
done_testing();
