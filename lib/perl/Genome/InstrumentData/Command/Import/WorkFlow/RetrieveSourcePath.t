#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

require Genome::Utility::Test;
require File::Temp;
use Test::More;

use_ok('Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePath') or die;
my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import') or die;

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $source_path = $test_dir.'/input.bam';
my $cmd = Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePath->execute(
    working_directory => $tmp_dir,
    source_path => $source_path,
);
ok($cmd, 'execute');
my $destination_path = $cmd->destination_path;
is($destination_path, $tmp_dir.'/input.bam', 'retrieved source path named correctly');
ok(-s $destination_path, 'destination path exists');

#print "$tmp_dir\n"; <STDIN>;
done_testing();
