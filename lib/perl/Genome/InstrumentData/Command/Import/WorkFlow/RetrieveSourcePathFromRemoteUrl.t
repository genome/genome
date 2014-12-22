#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require Genome::Utility::Test;
use File::Spec;
use File::Temp;
use Test::More;

use_ok('Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathFromRemoteUrl') or die;
my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'bam/v5') or die;
my $source_basename = 'input.bam';
my $source_path = 'file:/'.$test_dir.'/'.$source_basename;
my $source_md5_path = $source_path.'.md5';

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $cmd = Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathFromRemoteUrl->execute(
    working_directory => $tmp_dir,
    source_path => $source_path,
);
ok($cmd->result, 'execute');
my $destination_path = $cmd->destination_path;
is($destination_path, File::Spec->join($tmp_dir, $source_basename), 'destination path named correctly');
ok(-s $destination_path, 'destination path exists');

my $destination_md5_path = $cmd->destination_md5_path;
my $original_md5_path = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->original_md5_path_for($destination_path);
is($destination_md5_path, $original_md5_path, 'correcly named destination md5 path');
ok(-s $destination_md5_path, 'destination md5 path exists');
my $md5 = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->load_md5($destination_md5_path);
is($md5, 'f81fbc3d3a6b57d11e60b016bb2c950c', 'correct md5');

done_testing();
