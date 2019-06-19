#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use File::Spec;
use Test::More tests => 23;

my $class = 'Genome::Model::CwlPipeline::Command::Run';

use_ok($class);

my $build = Genome::Model::Build::CwlPipeline->__define__();
$build->data_directory( Genome::Sys->create_temp_directory );

my $cmd = $class->create(build => $build);

my ($workflow_dir, $results_dir) = $cmd->prepare_directories;
ok(-d $workflow_dir, 'created scratch directory');
ok(-d $results_dir, 'created results directory');

my @files;
my @paths;

for my $dir ('first', 'second', 'third') {
    Genome::Sys->create_directory(File::Spec->join($workflow_dir, $dir));
    for my $file (map { "${dir}${_}.txt" } (1..3)) {
        my $path = File::Spec->join($workflow_dir, $dir, $file);
        Genome::Sys->write_file($path, "hello $file\n");
        push @files, $file;
        push @paths, $path;
    }
}

my @test_structures = (
    [
        [
            { location => $paths[0], secondaryFiles => [ map { { location => $_ } } @paths[1,2] ] },
            { location => $paths[3], secondaryFiles => [] },
        ],
        [
            { location => $paths[4], secondaryFiles => [ { location => $paths[5] } ] },
        ],
    ],
    [
        { location => $paths[6], secondaryFiles => [] },
    ],
    { location => $paths[7], secondaryFiles => [ { location => $paths[8] } ] },
);

my $desired_prefix = 'turkey';
Genome::Model::Build::Input->create(
    build_id => $build->id,
    name => 'output_prefix',
    value_id => $desired_prefix,
    value_class_name => 'UR::Value::Text',
);
Genome::Model::Build::Input->create(
    build_id => $build->id,
    name => 'not_output_prefix',
    value_id => 'ignored',
    value_class_name => 'UR::Value::Text',
);


my $prefix = $cmd->_determine_output_prefix;
is($prefix, $desired_prefix, 'determined prefix');

for my $test (@test_structures) {
    $cmd->_stage_cromwell_output($results_dir, $test, $prefix);
}

for my $path (@paths) {
    ok(!-e $path, 'file was moved out');
}
for my $file (@files) {
    my $dest = File::Spec->join($results_dir, "$prefix-$file");
    ok(-e $dest, 'file was moved in with correct prefix applied');
}

$cmd->cleanup;
ok(!-d $workflow_dir, 'cleaned up scratch directory');

