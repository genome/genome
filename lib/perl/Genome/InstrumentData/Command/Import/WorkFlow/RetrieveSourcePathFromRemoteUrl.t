#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require Genome::Utility::Test;
use Test::More;

use_ok('Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathFromRemoteUrl') or die;
my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'bam/v5') or die;
my $source_path = 'file://'.$test_dir.'/input.bam';
my $source_md5_path = $test_dir.'/input.bam.md5';

my $cmd = Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathFromRemoteUrl->execute(
    working_directory => $test_dir,
    source_path => $source_path,
);
ok($cmd->result, 'execute');
my $destination_path = $cmd->destination_path;
is($destination_path, $test_dir.'/input.bam', 'destination path named correctly');
ok(-s $destination_path, 'destination path exists');
my $original_md5_path = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->original_md5_path_for($destination_path);
ok(-s $original_md5_path, 'destination md5 exists');

done_testing();
