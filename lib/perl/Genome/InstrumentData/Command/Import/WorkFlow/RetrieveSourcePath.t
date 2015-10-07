#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require Genome::Utility::Test;
require File::Spec;
use Test::More;

use_ok('Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePath') or die;
my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'v1') or die;

# Minimal base class testing
class Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathTest {
    is => 'Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePath',
};
sub Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathTest::execute { return 1; };

my $tmpdir = File::Spec->join(File::Spec->rootdir, 'tmp');
my $cmd = Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathTest->execute(
    working_directory => $tmpdir,
    source_path => 'input.bam',
);
ok($cmd->result, 'execute');
my $destination_file = Genome::InstrumentData::Command::Import::WorkFlow::SourceFile->create(path => File::Spec->join($tmpdir, 'input.bam'));
is($cmd->destination_path, $destination_file->path, 'destination path named correctly');
is($cmd->destination_original_md5_path, $destination_file->original_md5_path, 'destination md5 path named correctly');

done_testing();
