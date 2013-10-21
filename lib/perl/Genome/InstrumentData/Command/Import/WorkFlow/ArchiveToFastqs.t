#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

require Genome::Utility::Test;
require File::Compare;
require File::Temp;
use Test::More;

use_ok('Genome::InstrumentData::Command::Import::WorkFlow::ArchiveToFastqs') or die;
my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'v1') or die;
my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);

my $archive_base_name = 'input.2fastq.tgz';
my $archive_path = $tmp_dir.'/'.$archive_base_name;
Genome::Sys->create_symlink($test_dir.'/'.$archive_base_name, $archive_path);

my $cmd = Genome::InstrumentData::Command::Import::WorkFlow::ArchiveToFastqs->execute(
    working_directory => $tmp_dir,
    archive_path => $archive_path,
);
ok($cmd->result, 'execute');
my @fastq_paths = $cmd->fastq_paths;
for (my $i = 0; $i < @fastq_paths; $i++) {
    is($fastq_paths[$i], $tmp_dir.'/input.'.($i+1).'.fastq', "fastq path named correctly");
    ok(-s $fastq_paths[$i], 'fastq path exists');
    is(File::Compare::compare($fastq_paths[$i], $test_dir.'/input.'.($i + 1).'.fastq'), 0, 'fastq matches');
}

#print "$tmp_dir\n"; <STDIN>;
done_testing();
