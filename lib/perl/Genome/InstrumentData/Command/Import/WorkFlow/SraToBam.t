#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_COMMAND_DUMP_DEBUG_MESSAGES} = 1;
    $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
}

use strict;
use warnings;

use above "Genome";

require File::Compare;
require File::Temp;
require Genome::Utility::Test;
require Path::Class;
use Test::More;

use_ok('Genome::InstrumentData::Command::Import::Basic') or die;

my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'sra/v4');
my $sra_basename = 'input.sra';
my $source_sra_path = Path::Class::File->new($test_dir, $sra_basename);
ok(-s $source_sra_path, 'source sra exists') or die;
my $expected_bam_path = Path::Class::File->new($test_dir, $sra_basename.'.bam');
ok(-s $expected_bam_path, 'expected bam exists') or die;

my $tempdir = File::Temp::tempdir(CLEANUP => 1);
my $sra_path = Path::Class::File->new($tempdir, $sra_basename);
Genome::Sys->create_symlink($source_sra_path, $sra_path);

my $sample_name = '__TEST_SAMPLE__';
my $library_name = join('-', $sample_name, 'extlibs');

my $cmd = Genome::InstrumentData::Command::Import::WorkFlow::SraToBam->execute(
    sra_path => "$sra_path", # stringify
    working_directory => $tempdir,
    sample_name => $sample_name,
    library_name => $library_name,
);
ok($cmd->result, "execute sra to bam");

my $bam_path = $cmd->output_bam_path;
ok(-s $bam_path, 'bam path exists');
is($bam_path, $sra_path.'.bam', 'bam path correctly named');
is(File::Compare::compare($bam_path.'.flagstat', $expected_bam_path.'.flagstat'), 0, 'flagstats match');

done_testing();
