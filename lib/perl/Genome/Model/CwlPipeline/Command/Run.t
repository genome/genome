#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use File::Spec;
use Test::More tests => 19;

my $class = 'Genome::Model::CwlPipeline::Command::Run';

use_ok($class);

my $workflow_dir = Genome::Sys->create_temp_directory();
my $results_dir = Genome::Sys->create_temp_directory();


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


my $build = Genome::Model::Build::CwlPipeline->__define__();

my $cmd = $class->create(build => $build);

for my $test (@test_structures) {
    $cmd->_stage_cromwell_output($results_dir, $test);
}

for my $path (@paths) {
    ok(!-e $path, 'file was moved out');
}
for my $file (@files) {
    my $dest = File::Spec->join($results_dir, $file);
    ok(-e $dest, 'file was moved in');
}
