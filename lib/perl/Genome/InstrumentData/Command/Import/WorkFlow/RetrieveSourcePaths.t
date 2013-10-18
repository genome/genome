#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

require Genome::Utility::Test;
require File::Temp;
use Test::More;

use_ok('Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePath') or die;
my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'v1') or die;

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $source_path = $test_dir.'/input.bam';
my $cmd = Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePaths->execute(
    working_directory => $tmp_dir,
    source_paths => [$source_path],
);
ok($cmd, 'execute');
my @destination_paths = $cmd->destination_paths;
is_deeply(\@destination_paths, [$tmp_dir.'/input.bam'], 'destination paths named correctly');
ok(-s $destination_paths[0], 'destination path exists');
ok(-s $destination_paths[0].'.md5-orig', 'destination md5 exists');

#print "$tmp_dir\n"; <STDIN>;
done_testing();
