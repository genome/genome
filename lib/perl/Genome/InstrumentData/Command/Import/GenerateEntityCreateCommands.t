#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use strict;
use warnings;

use above 'Genome';

require File::Spec;
use Genome::Utility::Test;
use Test::Exception;
use Test::More;

use_ok('Genome::InstrumentData::Command::Import::GenerateEntityCreateCommands') or die;
my $data_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'generate-cmds');
my $input_file = File::Spec->join($data_dir, 'info.csv');

my $tmpdir = Genome::Sys->create_temp_directory;
my $commands_file_name = 'entity-commands-output.1.sh';
my $commands_file = File::Spec->join($tmpdir, $commands_file_name);
my $cmd = Genome::InstrumentData::Command::Import::GenerateEntityCreateCommands->execute(
    file => $input_file,
    output_file => $commands_file,
);
ok($cmd->result, 'execute');
Genome::Utility::Test::compare_ok($cmd->output_file, File::Spec->join($data_dir, $commands_file_name), 'commands files matches');

# Define Bwambale [individual create command should not be present]
my $bwambale = Genome::Individual->__define__(name => 'TGI-Bwambale');
ok($bwambale, '__define__ bwambale');

$commands_file_name = 'entity-commands-output.2.sh';
$commands_file = File::Spec->join($tmpdir, $commands_file_name);
$cmd = Genome::InstrumentData::Command::Import::GenerateEntityCreateCommands->execute(
    file => $input_file,
    output_file => $commands_file,
);
ok($cmd->result, 'execute');
Genome::Utility::Test::compare_ok($cmd->output_file, File::Spec->join($data_dir, $commands_file_name), 'commands files matches');

# Define Bwambale sample and library [all create command should not be present]
my $bwambale_sample = Genome::Sample->__define__(name => $bwambale->name.'-SRS394801', source => $bwambale);
ok($bwambale_sample, '__define__ bwambale sample');
my $bwambale_library = Genome::Library->__define__(name => $bwambale_sample->name.'-extlibs', sample => $bwambale_sample);
ok($bwambale_library, '__define__ bwambale library');

$commands_file_name = 'entity-commands-output.3.sh';
$commands_file = File::Spec->join($tmpdir, $commands_file_name);
$cmd = Genome::InstrumentData::Command::Import::GenerateEntityCreateCommands->execute(
    file => $input_file,
    output_file => $commands_file,
);
ok($cmd->result, 'execute');
Genome::Utility::Test::compare_ok($cmd->output_file, File::Spec->join($data_dir, $commands_file_name), 'commands files match');

# Fails
## missing required property
throws_ok(
    sub{ Genome::InstrumentData::Command::Import::GenerateEntityCreateCommands->execute(file => File::Spec->join($data_dir, 'missing-required-property.csv')); },
    qr/Missing required individual properties\: taxon/,
    'failed when missing requried property [individual taxon]',
);

## unknown property
throws_ok(
    sub{ Genome::InstrumentData::Command::Import::GenerateEntityCreateCommands->execute(file => File::Spec->join($data_dir, 'unknown-property.csv')); },
    qr/Unknown individual property\: unknown_property/,
    'failed w/ unknown property',
);

done_testing();
