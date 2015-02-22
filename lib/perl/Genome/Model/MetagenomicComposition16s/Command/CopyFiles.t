#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Genome::Model::Build::MetagenomicComposition16s::TestBuildFactory;
use Test::More;
use File::Temp;

BEGIN {
    use_ok('Genome::Model::MetagenomicComposition16s::Command::CopyFiles');
};

my ($build, $example_build);
subtest 'setup test builds' => sub {
    ($build, $example_build) = Genome::Model::Build::MetagenomicComposition16s::TestBuildFactory->build_with_example_build_for_sanger;
    ok($example_build, 'example build') or die;
    ok($example_build->get_or_create_data_directory, 'resolved data dir');

    is(        $build->status('Succeeded'), 'Succeeded', 'build is succeeded');
    is($example_build->status('Succeeded'), 'Succeeded', 'example_build is succeeded');

    my $time = time();
    my $date_template = UR::Context->date_template;
    my $older_date    = Date::Format::time2str($date_template, $time - 60);
    my $newer_date    = Date::Format::time2str($date_template, $time + 60);
    ok(        $build->date_completed($older_date), 'set build date_completed');
    ok($example_build->date_completed($newer_date), 'set example_build date_completed');
};

my $model = $example_build->model;
ok($model, 'got model from example_build');

# We run the "same" command three times so $tmdir and %copy_files_params are
# shared with those three subtests.
my $tmpdir = tmpdir();
my %copy_files_params = (
    models => [$model],
    builds => [$example_build],
    file_type => 'processed_fasta',
    destination => $tmpdir,
);

subtest 'run command with models and builds' => sub {
    my $cmd = Genome::Model::MetagenomicComposition16s::Command::CopyFiles->create(%copy_files_params);
    ok($cmd, 'created command with common params');
    $cmd->dump_status_messages(1);
    ok($cmd->execute, 'executed command');

    my @files = glob("$tmpdir/*");
    is(scalar @files, 1, 'Copied files');
};

subtest 're-run command to verify it fails due to existing output' => sub {
    my $cmd = Genome::Model::MetagenomicComposition16s::Command::CopyFiles->create(%copy_files_params);
    ok($cmd, 'created command with common params');
    $cmd->dump_status_messages(1);
    ok(!$cmd->execute, 'execute failed as expected due to existing output');
};

subtest 're-run command with --force to verify it overwrites existing output' => sub {
    my $cmd = Genome::Model::MetagenomicComposition16s::Command::CopyFiles->create(
        %copy_files_params,
        force => 1,
    );
    ok($cmd, 'created command with common params');
    $cmd->dump_status_messages(1);
    ok($cmd->execute, 'executed with --force option');
};

subtest 'run command with duplicate models and no builds' => sub {
    my $tmpdir = tmpdir();
    my $cmd = Genome::Model::MetagenomicComposition16s::Command::CopyFiles->create(
        models => [$model, $model],
        file_type => 'processed_fasta',
        destination => $tmpdir,
    );
    ok($cmd, 'created command');
    $cmd->dump_status_messages(1);
    ok($cmd->execute, 'executed command');
};

subtest 'run command with builds instead of models' => sub {
    my $tmpdir = tmpdir();
    my $cmd = Genome::Model::MetagenomicComposition16s::Command::CopyFiles->create(
        builds => [$example_build],
        file_type => 'processed_fasta',
        destination => $tmpdir,
    );
    ok($cmd, 'created command');
    $cmd->dump_status_messages(1);
    ok($cmd->execute, 'executed command');
};

subtest 'run command with missing file_type to verify failure' => sub {
    my $cmd = Genome::Model::MetagenomicComposition16s::Command::CopyFiles->create(
        builds => [$example_build],
    );
    ok($cmd, 'created command without file_type');
    $cmd->dump_status_messages(1);
    ok(!$cmd->execute, 'execute failed without file_type');
};

subtest 'run command with invalid file_type to verify failure' => sub {
    my $cmd = Genome::Model::MetagenomicComposition16s::Command::CopyFiles->create(
        models => [$model],
        file_type => 'some_file_type_that_is_not_valid',
    );
    ok($cmd, 'created command with invalid file_type');
    $cmd->dump_status_messages(1);
    ok(!$cmd->execute, 'execute failed with invalid file_type');
};

done_testing();

sub tmpdir {
    my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
    unless ($tmpdir && -d $tmpdir) {
        die 'failed to create tmpdir';
    }
    return $tmpdir;
}
