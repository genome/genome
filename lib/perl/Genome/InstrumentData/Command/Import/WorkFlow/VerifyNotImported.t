#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require Genome::InstrumentData::Command::Import::Inputs::SourceFile;
require Genome::Utility::Test;
require File::Spec;
require File::Temp;
use Test::More;

use_ok('Genome::InstrumentData::Command::Import::WorkFlow::VerifyNotImported') or die;
my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'bam/v1') or die;
my $source_file_basename = 'input.bam';
my $tempdir1 = File::Temp::tempdir(CLEANUP => 1);
my $tempdir2 = File::Temp::tempdir(CLEANUP => 1);

# Run MD5
my $source_path = File::Spec->join($tempdir1, $source_file_basename);
symlink(File::Spec->join($test_dir, $source_file_basename),  $source_path);
my $cmd = Genome::InstrumentData::Command::Import::WorkFlow::VerifyNotImported->create(
    working_directory => $tempdir1,
    source_paths => [$source_path],
);
is_deeply([$cmd->output_path], [$cmd->output_paths], 'output_path and output_paths');
my $source_file = $cmd->source_file_for_path($source_path);
is($source_file->path, $source_path, 'source_file');
my $output_file = $cmd->output_file_for_path($source_path);
is($output_file->path, File::Spec->join($tempdir1, File::Basename::basename($source_path)), 'output_file');
ok(!$output_file->md5_path_size, 'MD5 path does not exist');
ok($cmd->execute, 'execute');
ok($output_file->md5_path_size, 'MD5 path exists');

# Copy the MD5 that was created above
$cmd = Genome::InstrumentData::Command::Import::WorkFlow::VerifyNotImported->create(
    working_directory => $tempdir2,
    source_paths => [$source_path],
);
$output_file = $cmd->output_file_for_path($source_path);
ok(!$output_file->md5_path_size, 'MD5 path does not exist');
ok($cmd->execute, 'execute');
ok($output_file->md5_path_size, 'MD5 path exists');
unlink $output_file->md5_path; # remove since we are using the same directory below

# Previously Imported MD5
my $instdata = Genome::InstrumentData::Imported->__define__(id => -11);
my $md5_attr = Genome::InstrumentDataAttribute->__define__(
    instrument_data_id => $instdata->id,
    attribute_label => 'original_data_path_md5',
    attribute_value => '940825168285c254b58c47399a3e1173',
    nomenclature => 'WUGC',
);
ok($md5_attr, 'create md5 inst data attr');
$cmd = Genome::InstrumentData::Command::Import::WorkFlow::VerifyNotImported->execute(
    working_directory => $tempdir2,
    source_paths => [$source_path],
);
ok(!$cmd->result, 'execute');
is(Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get->error_message, 'Instrument data was previously imported! Found existing instrument data: -11', 'correct error');
$output_file = $cmd->output_file_for_path($source_path);
unlink $output_file->md5_path;

# With downsampling
$cmd = Genome::InstrumentData::Command::Import::WorkFlow::VerifyNotImported->execute(
    working_directory => $tempdir2,
    source_paths => [$source_path],
    downsample_ratio => 0.25,
);
ok($cmd->result, 'execute');
$output_file = $cmd->output_file_for_path($source_path);
unlink $output_file->md5_path;

# previously imported...
my $downsample_ratio_attr = Genome::InstrumentDataAttribute->__define__(
    instrument_data_id => -11,
    attribute_label => 'downsample_ratio',
    attribute_value => 0.25,
    nomenclature => 'WUGC',
);
ok($downsample_ratio_attr, '__define__ downsample_ratio_attr');
$cmd = Genome::InstrumentData::Command::Import::WorkFlow::VerifyNotImported->execute(
    working_directory => $tempdir2,
    source_paths => [$source_path],
    downsample_ratio => .25,
);
ok(!$cmd->result, 'execute');
is(Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get->error_message, 'Instrument data was previously downsampled by a ratio of 0.25 and imported! Found existing instrument data: -11', 'correct error');

#print "$tmp_dir\n"; <STDIN>;
done_testing();
