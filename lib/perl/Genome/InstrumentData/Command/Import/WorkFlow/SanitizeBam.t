#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

require Genome::Utility::Test;
require File::Temp;
require File::Compare;
use Test::More;

use_ok('Genome::InstrumentData::Command::Import::WorkFlow::SanitizeBam') or die;
my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'bam/v1') or die;

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $dirty_bam_path = $tmp_dir.'/input.bam';
Genome::Sys->create_symlink($test_dir.'/input.bam', $dirty_bam_path);
ok(-s $dirty_bam_path, 'linked dirty bam path');
Genome::Sys->create_symlink($test_dir.'/input.bam.md5', $dirty_bam_path.'.md5');
ok(-s $dirty_bam_path.'.md5', 'linked dirty bam md5 path');
Genome::Sys->create_symlink($test_dir.'/input.bam.flagstat', $dirty_bam_path.'.flagstat');
ok(-s $dirty_bam_path.'.flagstat', 'linked dirty bam flagstat path');

my $cmd = Genome::InstrumentData::Command::Import::WorkFlow::SanitizeBam->execute(
    dirty_bam_path => $dirty_bam_path,
);
ok($cmd, 'execute');
my $clean_bam_path = $cmd->clean_bam_path;
is($clean_bam_path, $tmp_dir.'/input.clean.bam', 'clean bam path named correctly');
ok(-s $clean_bam_path, 'clean bam path exists');
is(File::Compare::compare($clean_bam_path, $test_dir.'/input.clean.bam'), 0, 'clean bam matches');
ok(-s $clean_bam_path.'.flagstat', 'flagstat path exists');
is(File::Compare::compare($clean_bam_path.'.flagstat', $test_dir.'/input.clean.bam.flagstat'), 0, 'flagstat matches');

ok(!-e $dirty_bam_path, 'removed dirty bam path after sanitizing');
ok(!-e $dirty_bam_path.'.md5', 'removed dirty md5 path after sanitizing');
ok(!-e $dirty_bam_path.'.md5-orig', 'removed dirty md5 orig path after sanitizing');
ok(!-e $dirty_bam_path.'.flagstat', 'removed dirty flagstat path after sanitizing');

#print "$tmp_dir\n"; <STDIN>;
done_testing();
