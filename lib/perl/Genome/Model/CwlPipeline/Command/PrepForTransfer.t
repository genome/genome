#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
      $ENV{UR_DBI_NO_COMMIT} = 1;
}

use above 'Genome';
use Genome::Test::Factory::Sample;
use Genome::Test::Factory::Build;
use Genome::Test::Factory::Model::CwlPipeline;
use Genome::Utility::Test 'compare_ok';
use Test::More;
use File::Spec qw();



my $pkg = 'Genome::Model::CwlPipeline::Command::PrepForTransfer';
use_ok($pkg);

my $data_dir = sprintf "%s.d", __FILE__;

my @dirs = qw/A B C D/;
my @builds;
for my $dir (@dirs) {
    my $sample = Genome::Test::Factory::Sample->setup_object(
       name => 'sample-'. $dir,
    );
    my $model = Genome::Test::Factory::Model::CwlPipeline->setup_object(
       name => 'test-model-'. $dir,
       subject => $sample,
    );
    my $build = Genome::Test::Factory::Build->setup_object(
       id => 'build'. $dir,
       model_id => $model->id,
       data_directory => File::Spec->catfile($data_dir,$dir),
    );
    push @builds, $build;
}

my $expected_file = File::Spec->catfile($data_dir, 'MANIFEST');

my $output_dir = Genome::Sys->create_temp_directory;
my $output_file = File::Spec->catfile($output_dir,'MANIFEST');

my $cmd = $pkg->create(
   builds => \@builds,
   directory => $output_dir,
   md5sum => 1,
);

ok($cmd, "created command");
ok($cmd->execute, "executed command");

compare_ok($expected_file, $output_file);

done_testing();
