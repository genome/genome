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

use_ok('Genome::InstrumentData::Command::Import::WorkFlow::DownsampleBam') or die;

my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'bam/v3') or die;

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $bam_base_name = 'test.bam';
my $bam_path = $tmp_dir.'/'.$bam_base_name;
Genome::Sys->create_symlink($test_dir.'/'.$bam_base_name, $bam_path);
ok(-s $bam_path, 'linked bam path');
Genome::Sys->create_symlink($test_dir.'/'.$bam_base_name.'.flagstat', $bam_path.'.flagstat');
ok(-s $bam_path.'.flagstat', 'linked bam flagstat path');

# Check downsample value errors
my $cmd = Genome::InstrumentData::Command::Import::WorkFlow::DownsampleBam->create(
    bam_path => $bam_path,
    downsample_ratio => 'NA',
);
ok($cmd, 'create');
my @errors = $cmd->__errors__;
ok(@errors, 'errors for downsample_ratio of NA');
is($errors[0]->desc, 'Invalid number! NA', 'correct error desc for downsample_ratio of NA');

$cmd = Genome::InstrumentData::Command::Import::WorkFlow::DownsampleBam->create(
    bam_path => $bam_path,
    downsample_ratio => 0,
);
ok($cmd, 'create');
@errors = $cmd->__errors__;
ok(@errors, 'errors for downsample_ratio of 0');
is($errors[0]->desc, 'Must be greater than 0 and less than 1! 0', 'correct error desc for downsample_ratio of 0');

$cmd = Genome::InstrumentData::Command::Import::WorkFlow::DownsampleBam->create(
    bam_path => $bam_path,
    downsample_ratio => 1,
);
ok($cmd, 'create');
@errors = $cmd->__errors__;
ok(@errors, 'errors for downsample_ratio of 1');
is($errors[0]->desc, 'Must be greater than 0 and less than 1! 1', 'correct error desc for downsample_ratio of 1');

# Success
$cmd = Genome::InstrumentData::Command::Import::WorkFlow::DownsampleBam->create(
    bam_path => $bam_path,
    downsample_ratio => .25,
);
ok($cmd, 'create');
ok($cmd->execute, 'execute');
my $downsampled_bam_path = $cmd->downsampled_bam_path;
my $downsampled_bam_base_name = 'test.downsampled.bam';
is($downsampled_bam_path, $tmp_dir.'/'.$downsampled_bam_base_name, 'downsampled bam path named correctly');
ok(-s $downsampled_bam_path, 'downsampled bam path exists');
is(File::Compare::compare($downsampled_bam_path, $test_dir.'/'.$downsampled_bam_base_name), 0, 'downsampled bam matches');
ok(-s $downsampled_bam_path.'.flagstat', 'flagstat path exists');
is(File::Compare::compare($downsampled_bam_path.'.flagstat', $test_dir.'/'.$downsampled_bam_base_name.'.flagstat'), 0, 'flagstat matches');
ok(!glob($bam_path.'*'), 'removed bam path and auxillary files after down sampling');

#print "$tmp_dir\n"; <STDIN>;
done_testing();
