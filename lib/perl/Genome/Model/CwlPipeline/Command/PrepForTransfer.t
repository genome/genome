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
my %files_by_build;
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
    $build->date_completed('1999-04-0' . (ord($dir)-60) . ' 18:04:01');
    push @builds, $build;
    my @files = grep {-f $_} glob($build->data_directory .'/results/*');
    $files_by_build{$build->id} = \@files;
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

my $rv;
eval {
   $rv = $cmd->execute;
};    

ok($rv, "executed command");

compare_ok($expected_file, $output_file);

for my $build_id (sort keys %files_by_build) {
    for my $file (@{$files_by_build{$build_id}}) {
        my ($file_name, $dir, $suffix) = File::Basename::fileparse($file);
        my $symlink_name = $build_id .'.'. $file_name;
        my $symlink_path = File::Spec->join($output_dir,$symlink_name); 
        ok(-l $symlink_path, 'symlink exists');
        my $target = readlink($symlink_path);
        ok(-e $target, 'target exists');        
    }
}   


my $cmd2 = $pkg->create(
   builds => \@builds,
   directory => $output_dir,
   md5sum => 1,
);

ok($cmd2, "created command");

my $rv2;
eval {
     $rv2 = $cmd2->execute;
};    
ok(!$rv2, "failed to execute command a second time");

done_testing();
