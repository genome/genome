#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require Genome::Utility::Test;
require File::Spec;
require File::Temp;
use Test::More;

use_ok('Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathFromLocalDisk') or die;
my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'bam/v5') or die;
my $source_basename = 'input.bam';
my $source_path =File::Spec->join($test_dir, $source_basename);
my $source_md5_path = $source_path.'.md5';

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $cmd = Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathFromLocalDisk->execute(
    working_directory => $tmp_dir,
    source_path => $source_path,
);
ok($cmd->result, 'execute');
my $destination_path = $cmd->destination_path;
is($destination_path, File::Spec->join($tmp_dir, $source_basename), 'destination path named correctly');
ok(-s $destination_path, 'destination path exists');

ok(-s $cmd->destination_path, 'destination path exists');
like($cmd->destination_original_md5_path, qr/\.md5-orig$/, 'correctly named destination_original_md5_path');
ok(-s $cmd->destination_original_md5_path, 'destination_original_md5_path exists');
my $md5 = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->load_md5($cmd->destination_original_md5_path);
is($md5, 'f81fbc3d3a6b57d11e60b016bb2c950c', 'correct md5');

done_testing();
