#! /gsc/bin/perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
}

use strict;
use warnings;

use above 'Genome';

require Genome::Utility::Test;
require File::Compare;
require File::Spec;
require File::Temp;
use Test::More;

use_ok('Genome::InstrumentData::Command::Import::WorkFlow::DownsampleBam') or die;

my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', File::Spec->catfile('bam', 'v4')) or die;
my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $bam_base_name = 'test.clean.sorted.bam';
my $bam_path = File::Spec->catfile($tmp_dir, $bam_base_name);
Genome::Sys->create_symlink( File::Spec->catfile($test_dir, $bam_base_name), $bam_path);
ok(-s $bam_path, 'linked bam path');
Genome::Sys->create_symlink( File::Spec->catfile($test_dir, $bam_base_name.'.flagstat'), $bam_path.'.flagstat');
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
my $output_bam_path = $cmd->output_bam_path;
my $downsampled_bam_base_name = 'test.clean.sorted.downsampled.bam';
is($output_bam_path, File::Spec->catfile($tmp_dir, $downsampled_bam_base_name), 'downsampled bam path named correctly');
ok(-s $output_bam_path, 'downsampled bam path exists');
is(File::Compare::compare($output_bam_path, File::Spec->catfile($test_dir, $downsampled_bam_base_name)), 0, 'downsampled bam matches');
ok(-s $output_bam_path.'.flagstat', 'flagstat path exists');
is(File::Compare::compare($output_bam_path.'.flagstat', File::Spec->catfile($test_dir, $downsampled_bam_base_name.'.flagstat')), 0, 'flagstat matches');
ok(!glob($bam_path.'*'), 'removed bam path and auxillary files after down sampling');

#print "$tmp_dir\n"; <STDIN>;
done_testing();
