#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';

use File::Spec;
use Test::More tests => 5;

my $class = 'Genome::ProcessingProfile::Command::Create::CwlPipeline';
use_ok($class);

my $dir = Genome::Sys->create_temp_directory;
my $filename = 'main.cwl';
my $file = File::Spec->join($dir, $filename);

Genome::Sys->write_file($file, 'fake');

my $cmd = $class->create(
    name => 'test for cwl-pipeline pp create',
    cwl_directory => $dir,
    main_workflow_file => $filename,
    primary_docker_image => 'docker(mgibio/tabix)',
);
isa_ok($cmd, $class, 'created command');

ok($cmd->execute, 'executed command');
my $pp = $cmd->created_processing_profile;

isa_ok($pp, $class->_target_class_name, 'created processing profile');

my $diff = Genome::Sys->diff_file_vs_file($file, $pp->main_workflow_file);
ok(!$diff, 'file faithfully copied')
    or diag('diff: ' . $diff);
