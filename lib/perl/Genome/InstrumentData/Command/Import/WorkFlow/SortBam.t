#! /gsc/bin/perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
}

use strict;
use warnings;

use above 'Genome';

require Genome::Utility::Test;
require File::Temp;
require File::Compare;
use Test::More;

use_ok('Genome::InstrumentData::Command::Import::WorkFlow::SortBam') or die;
my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'bam/v1') or die;

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $unsorted_bam_base_name = 'input.clean.bam';
my $unsorted_bam_path = $tmp_dir.'/'.$unsorted_bam_base_name;
Genome::Sys->create_symlink($test_dir.'/'.$unsorted_bam_base_name, $unsorted_bam_path);
ok(-s $unsorted_bam_path, 'linked unsorted bam path');
#Genome::Sys->create_symlink($test_dir.'/'.$unsorted_bam_base_name.'.md5', $unsorted_bam_path.'.md5');
#ok(-s $unsorted_bam_path.'.md5', 'linked unsorted bam md5 path');
Genome::Sys->create_symlink($test_dir.'/'.$unsorted_bam_base_name.'.flagstat', $unsorted_bam_path.'.flagstat');
ok(-s $unsorted_bam_path.'.flagstat', 'linked unsorted bam flagstat path');

my $cmd = Genome::InstrumentData::Command::Import::WorkFlow::SortBam->execute(
    unsorted_bam_path => $unsorted_bam_path,
);
ok($cmd, 'execute');
my $sorted_bam_path = $cmd->sorted_bam_path;
my $sorted_bam_base_name = 'input.clean.sorted.bam';
is($sorted_bam_path, $tmp_dir.'/'.$sorted_bam_base_name, 'sorted bam path named correctly');
ok(-s $sorted_bam_path, 'sorted bam path exists');
is(File::Compare::compare($sorted_bam_path, $test_dir.'/'.$sorted_bam_base_name), 0, 'sorted bam matches');
ok(-s $sorted_bam_path.'.flagstat', 'flagstat path exists');
is(File::Compare::compare($sorted_bam_path.'.flagstat', $test_dir.'/'.$sorted_bam_base_name.'.flagstat'), 0, 'flagstat matches');

ok(!-e $unsorted_bam_path, 'removed unsorted bam path after sorting');
ok(!-e $unsorted_bam_path.'.md5', 'removed unsorted md5 path after sorting');
ok(!-e $unsorted_bam_path.'.md5-orig', 'removed unsorted md5 orig path after sorting');
ok(!-e $unsorted_bam_path.'.flagstat', 'removed unsorted flagstat path after sorting');

#print "$tmp_dir\n"; <STDIN>;
done_testing();
