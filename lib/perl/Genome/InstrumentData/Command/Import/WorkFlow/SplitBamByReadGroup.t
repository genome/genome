#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

require Genome::Utility::Test;
require File::Temp;
use Test::More;

use_ok('Genome::InstrumentData::Command::Import::WorkFlow::SplitBamByReadGroup') or die;
my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import') or die;

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $one_read_group_basename = 'input.bam';
my $one_read_group_bam_path = $tmp_dir.'/'.$one_read_group_basename;
Genome::Sys->create_symlink($test_dir.'/'.$one_read_group_basename, $one_read_group_bam_path);
ok(-s $one_read_group_bam_path, 'linked one read group bam');
my $cmd = Genome::InstrumentData::Command::Import::WorkFlow::SplitBamByReadGroup->execute(bam_path => $one_read_group_bam_path);
ok($cmd, 'execute');
my @read_group_bam_paths = $cmd->read_group_bam_paths;
is(@read_group_bam_paths, 1, 'got one read group bam path');
is($read_group_bam_paths[0], $tmp_dir.'/input.bam', 'copied single read group bam to tmp');
ok(-s $read_group_bam_paths[0], 'single read group bam exists');

my $multi_rg_base_name = 'input.rg-multi';
my $multi_rg_bam_path = $tmp_dir.'/'.$multi_rg_base_name.'.bam';
Genome::Sys->create_symlink($test_dir.'/'.$multi_rg_base_name.'.bam', $multi_rg_bam_path);
Genome::Sys->create_symlink($test_dir.'/'.$multi_rg_base_name.'bam..flagstat', $multi_rg_bam_path.'.flagstat');
ok(-s $multi_rg_bam_path, 'linked two read groups bam');
$cmd = Genome::InstrumentData::Command::Import::WorkFlow::SplitBamByReadGroup->execute(bam_path => $multi_rg_bam_path);
ok($cmd, 'execute');
@read_group_bam_paths = $cmd->read_group_bam_paths;
is(@read_group_bam_paths, 2, '2 read group bam paths');
is_deeply(
    \@read_group_bam_paths, 
    [ map { $tmp_dir.'/'.$multi_rg_base_name.'.'.$_.'.bam' } (qw/ 2883581797 2883581798 /)],
    '2 read groups bams paths');
is(grep({ -s } @read_group_bam_paths), 2, 'read group bam paths exist');

#print "$tmp_dir\n"; <STDIN>;
done_testing();
