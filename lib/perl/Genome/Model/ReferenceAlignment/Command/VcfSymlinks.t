#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Test::Factory::Model::ReferenceAlignment;
use Genome::Test::Factory::Build;
use IO::File;
use File::Spec;

# This test was auto-generated because './Model/ReferenceAlignment/Command/SubmissionSummary.pm'
# had no '.t' file beside it.  Please remove this test if you believe it was
# created unnecessarily.  This is a bare minimum test that just compiles Perl
# and the UR class.
use_ok('Genome::Model::ReferenceAlignment::Command::VcfSymlinks');

# how to test the rest of this? All it does is just create symlinks
my $model = Genome::Test::Factory::Model::ReferenceAlignment->setup_object();
ok($model->isa("Genome::Model::ReferenceAlignment"), "Generated a reference alignment model");
ok($model->processing_profile->indel_detection_strategy('samtools'), "Updated indel detection strategy");

my $test_dir = Genome::Sys->create_temp_directory;
isnt($test_dir, undef, "Created a temp directory: $test_dir");

my $build = Genome::Test::Factory::Build->setup_object(
    model_id => $model->id,
    data_directory => $test_dir,
    status => 'Succeeded',
);
ok($build, "Created a test build");

#create the variants subdir and some fake vcf files
my $variants_dir_path = File::Spec->join($test_dir, "variants");
my $mkdir_result = mkdir "$variants_dir_path";
ok($mkdir_result, "Created variants directory in build data directory");

my $snvs_path = File::Spec->join($variants_dir_path, "snvs.vcf.gz");
my $sfh = IO::File->new($snvs_path,"w");
ok($sfh, "Created an empty snvs.vcf.gz");
$sfh->close;

my $indels_path = File::Spec->join($variants_dir_path, "indels.vcf.gz");
my $ifh = IO::File->new($indels_path,"w");
ok($ifh, "Created an empty snvs.vcf.gz");
$ifh->close;

my $output_test_dir = Genome::Sys->create_temp_directory;
isnt($output_test_dir, undef, "Created an output temp directory");

my $command = Genome::Model::ReferenceAlignment::Command::VcfSymlinks->execute(
    output_directory => $output_test_dir,
    models => [$model],
);

ok($command->result,"Successfully executed command");
my $sample = $model->subject_name;
isnt($sample, undef, "Sample name $sample retrieved from model");
ok(-e File::Spec->join($output_test_dir,"$sample.snvs.vcf.gz"), "SNVs VCF linked in");
ok(-e File::Spec->join($output_test_dir,"$sample.indels.vcf.gz"), "Indels VCF linked in");

done_testing();
