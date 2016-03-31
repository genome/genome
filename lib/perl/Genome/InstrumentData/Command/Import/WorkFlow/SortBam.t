#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
    $ENV{UR_COMMAND_DUMP_DEBUG_MESSAGES} = 1;
}

use strict;
use warnings;

use above 'Genome';

require Genome::Utility::Test;
require File::Temp;
require File::Compare;
use Test::More;

use_ok('Genome::InstrumentData::Command::Import::WorkFlow::SortBam') or die;
my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'v01') or die;

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $unsorted_bam_base_name = 'sort-bam.input.bam';
my $bam_path = $tmp_dir.'/'.$unsorted_bam_base_name;
Genome::Sys->create_symlink($test_dir.'/'.$unsorted_bam_base_name, $bam_path);
ok(-s $bam_path, 'linked unsorted bam path');
Genome::Sys->create_symlink($test_dir.'/'.$unsorted_bam_base_name.'.flagstat', $bam_path.'.flagstat');
ok(-s $bam_path.'.flagstat', 'linked unsorted bam flagstat path');

my $cmd = Genome::InstrumentData::Command::Import::WorkFlow::SortBam->execute(
    working_directory => $tmp_dir,
    bam_path => $bam_path,
);
ok($cmd->result, 'execute');
my $output_bam_path = $cmd->output_bam_path;
is($output_bam_path, $tmp_dir.'/'. 'sort-bam.input.sorted.bam', 'sorted bam path named correctly');
ok(-s $output_bam_path, 'sorted bam path exists');
my $expected_bam_basename = 'sort-bam.expected.bam';
is(File::Compare::compare($output_bam_path, $test_dir.'/'.$expected_bam_basename), 0, 'sorted bam matches');
ok(-s $output_bam_path.'.flagstat', 'flagstat path exists');
is(File::Compare::compare($output_bam_path.'.flagstat', $test_dir.'/'.$expected_bam_basename.'.flagstat'), 0, 'flagstat matches');
ok(!glob($bam_path.'*'), 'removed bam path and auxillary files after sorting');

#print "$tmp_dir\n"; <STDIN>;
done_testing();
