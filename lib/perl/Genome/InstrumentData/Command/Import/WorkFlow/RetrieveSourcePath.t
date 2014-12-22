#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require Genome::Utility::Test;
require File::Spec;
require File::Temp;
use Test::More;

use_ok('Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePath') or die;
my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'bam/v1') or die;

class Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathTest {
    is => 'Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePath',
};
my $tmpdir = File::Spec->join(File::Spec->rootdir, 'tmp');
my $cmd = Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathTest->create(
    working_directory => $tmpdir,
    source_path => 'input.bam',
);
my $destination_path = $cmd->destination_path;
is($destination_path, File::Spec->join($tmpdir, 'input.bam'), 'destination path named correctly');
is(
    $cmd->destination_md5_path,
    Genome::InstrumentData::Command::Import::WorkFlow::Helpers->original_md5_path_for($destination_path),
    'destination md5 path named correctly',
);

done_testing();
