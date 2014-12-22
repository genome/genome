#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require Genome::Utility::Test;
require File::Compare;
require File::Temp;
use Test::More;

use_ok('Genome::InstrumentData::Command::Import::WorkFlow::ArchiveToFastqs') or die;
use_ok('Genome::InstrumentData::Command::Import::WorkFlow::Helpers') or die;
my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'fastq/v1') or die;
my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);

my $archive_base_name = 'input.fastq.tgz';
my $archive_source_path = $test_dir.'/'.$archive_base_name;
my $archive_path = $tmp_dir.'/'.$archive_base_name;
Genome::Sys->create_symlink($archive_source_path, $archive_path) or die;

my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;
Genome::Sys->create_symlink(
    $helpers->md5_path_for($archive_source_path),
    $helpers->md5_path_for($archive_path),
) or die;
Genome::Sys->create_symlink(
    $helpers->md5_path_for($archive_source_path),
    $helpers->original_md5_path_for($archive_path),
) or die;

my $cmd = Genome::InstrumentData::Command::Import::WorkFlow::ArchiveToFastqs->execute(
    working_directory => $tmp_dir,
    archive_path => $archive_path,
);
ok($cmd->result, 'execute');
my @fastq_paths = $cmd->fastq_paths;
for (my $i = 0; $i < @fastq_paths; $i++) {
    is($fastq_paths[$i], $tmp_dir.'/input.'.($i+1).'.fastq', "fastq path named correctly");
    ok(-s $fastq_paths[$i], 'fastq path exists');
    is(File::Compare::compare($fastq_paths[$i], $test_dir.'/input.'.($i + 1).'.fastq'), 0, 'fastq '.($i + 1).' matches');
}

ok(!glob($archive_path.'*'), 'removed archived source path after extracting');
ok(!-e $cmd->extract_directory, 'removed extract diectory after extracting');

# FAIL - no fastqs in archive
$cmd = Genome::InstrumentData::Command::Import::WorkFlow::ArchiveToFastqs->execute(
    working_directory => $tmp_dir,
    archive_path => $test_dir.'/input.fail.no-fastqs.tar',
);
ok(!$cmd->result, 'execute failed');
is($cmd->error_message, 'No fastqs found in archive!', 'correct error');

#print "$tmp_dir\n"; <STDIN>;
done_testing();
