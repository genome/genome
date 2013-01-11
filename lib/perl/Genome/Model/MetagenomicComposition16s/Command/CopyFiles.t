#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Genome::Model::Build::MetagenomicComposition16s::TestBuildFactory;
use Test::More;

use_ok('Genome::Model::MetagenomicComposition16s::Command::CopyFiles') or die;

my ($build, $example_build) = Genome::Model::Build::MetagenomicComposition16s::TestBuildFactory->build_with_example_build_for_sanger;
ok($example_build, 'example build') or die;
my $model = $example_build->model;
ok($example_build->get_or_create_data_directory, 'resolved data dir');

is(        $build->the_master_event->event_status('Succeeded'), 'Succeeded', 'build is succeeded');
is($example_build->the_master_event->event_status('Succeeded'), 'Succeeded', 'example_build is succeeded');

my $time = time();
my $older_date = Date::Format::time2str(q(%Y-%m-%d %H:%M:%S), $time - 60);
my $newer_date = Date::Format::time2str(q(%Y-%m-%d %H:%M:%S), $time + 60);
ok(        $build->the_master_event->date_completed($older_date), 'set build date_completed');
ok($example_build->the_master_event->date_completed($newer_date), 'set example_build date_completed');

my $tmpdir = File::Temp::tempdir(CLEANUP => 1);

# ok - copy w/  models and builds
my $cmd = Genome::Model::MetagenomicComposition16s::Command::CopyFiles->create(
    models => [$model],
    builds => [$example_build],
    file_type => 'processed_fasta',
    destination => $tmpdir,
);
ok($cmd, 'create');
$cmd->dump_status_messages(1);
ok($cmd->execute, 'execute');
my @files = glob("$tmpdir/*");
is(scalar @files, 1, 'Copied files');

# fail - copy to existing
$cmd = Genome::Model::MetagenomicComposition16s::Command::CopyFiles->create(
    models => [$model],
    file_type => 'processed_fasta',
    destination => $tmpdir,
);
ok($cmd, 'create');
$cmd->dump_status_messages(1);
ok(!$cmd->execute, 'execute failed as expected to copy existing file');

# ok - force copy
$cmd = Genome::Model::MetagenomicComposition16s::Command::CopyFiles->create(
    models => [$model, $model],
    file_type => 'processed_fasta',
    destination => $tmpdir,
    force => 1,
);
ok($cmd, 'create');
$cmd->dump_status_messages(1);
ok($cmd->execute, 'execute w/ force copy');

# fail - no type
$cmd = Genome::Model::MetagenomicComposition16s::Command::CopyFiles->create(
    builds => [$example_build],
);
ok($cmd, 'create');
$cmd->dump_status_messages(1);
ok(!$cmd->execute, 'execute failed w/o type');

# fail - invalid type
$cmd = Genome::Model::MetagenomicComposition16s::Command::CopyFiles->create(
    models => [$model],
    file_type => 'some_file_type_that_is_not_valid',
);
ok($cmd, 'create');
$cmd->dump_status_messages(1);
ok(!$cmd->execute, 'execute failed w/ invalid type');

done_testing();
