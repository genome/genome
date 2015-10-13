#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require Genome::Utility::Test;
require File::Compare;
require File::Spec;
require File::Temp;
use Test::More;

use_ok('Genome::InstrumentData::Command::Import::WorkFlow::ArchiveToFastqs') or die;
use_ok('Genome::InstrumentData::Command::Import::WorkFlow::SourceFile') or die;
my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'v1') or die;
my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);

my $source_base_name = 'input.fastq.tgz';
my $source_file = Genome::InstrumentData::Command::Import::WorkFlow::SourceFile->create(
    path => File::Spec->join($test_dir, $source_base_name),
);
my $local_source_file = Genome::InstrumentData::Command::Import::WorkFlow::SourceFile->create(
    path => File::Spec->join($tmp_dir, $source_base_name),
);
Genome::Sys->create_symlink($source_file->path, $local_source_file->path) or die;

Genome::Sys->create_symlink($source_file->md5_path, $local_source_file->md5_path) or die;
Genome::Sys->create_symlink($source_file->md5_path, $local_source_file->original_md5_path) or die;

my $cmd = Genome::InstrumentData::Command::Import::WorkFlow::ArchiveToFastqs->execute(
    working_directory => $tmp_dir,
    archive_path => $local_source_file->path,
);
ok($cmd->result, 'execute');
my @fastq_paths = $cmd->fastq_paths;
for (my $i = 0; $i < @fastq_paths; $i++) {
    is($fastq_paths[$i], File::Spec->join($tmp_dir, 'input.'.($i+1).'.fastq'), "fastq path named correctly");
    ok(-s $fastq_paths[$i], 'fastq path exists');
    is(File::Compare::compare($fastq_paths[$i], File::Spec->join($test_dir, 'input.'.($i + 1).'.fastq')), 0, 'fastq '.($i + 1).' matches');
}

ok(!glob($local_source_file->path.'*'), 'removed archived source path after extracting');
ok(!-e $cmd->extract_directory, 'removed extract directory after extracting');

# FAIL - no fastqs in archive
$cmd = Genome::InstrumentData::Command::Import::WorkFlow::ArchiveToFastqs->execute(
    working_directory => $tmp_dir,
    archive_path => File::Spec->join($test_dir, 'archive-to-fastq.no-fastqs.tar'),
);
ok(!$cmd->result, 'execute failed');
is($cmd->error_message, 'No fastqs found in archive!', 'correct error');

#print "$tmp_dir\n"; <STDIN>;
done_testing();
