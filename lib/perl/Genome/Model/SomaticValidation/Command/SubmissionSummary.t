#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Utility::Test;

use_ok("Genome::Test::Factory::Model::SomaticValidation");
use_ok("Genome::Test::Factory::Model::ReferenceSequence");
use_ok("Genome::Test::Factory::Build");
use_ok('Genome::Model::SomaticValidation::Command::SubmissionSummary');

#my $base_dir = Genome::Utility::Test->data_dir_ok("Genome::Model::SomaticValidation::Command::SubmissionSummary", "v1");
my $bam_dir = Genome::Sys->create_temp_directory;
my $ref_model = Genome::Test::Factory::Model::ReferenceSequence->setup_object;
my $ref_build = Genome::Test::Factory::Build->setup_object(model_id => $ref_model->id);
my $result = Genome::InstrumentData::AlignmentResult::Merged->__define__(id => "-12355", output_dir => $bam_dir);
my $bam_path = $bam_dir."/".$result->id.".bam";
`touch $bam_path`;
my $validation_model = Genome::Test::Factory::Model::SomaticValidation->setup_object(reference_sequence_build => $ref_build, tumor_sample => Genome::Sample->create(name => "tumor1"));
my $validation_build = Genome::Test::Factory::Build->setup_object(model_id => $validation_model->id, status => "Succeeded");
my $u = Genome::SoftwareResult::User->create(label => "merged_alignment", software_result => $result, user => $validation_build);

my $output_dir = Genome::Sys->create_temp_directory;
my $samples_file = "$output_dir/sample_mapping";
my $bams_file = "$output_dir/bam_list";
my $md5s_file = "$output_dir/md5s";
my $command = Genome::Model::SomaticValidation::Command::SubmissionSummary->create(sample_mapping_file => $samples_file, bam_list_file => $bams_file, md5_list_file => $md5s_file, models => [$validation_model]);

ok($command->execute, "Executed successfully");
ok(-s $samples_file, "Samples file was created");
ok(-s $bams_file, "Bams file was created");
ok(-s $md5s_file, "md5 file was created");

done_testing;
