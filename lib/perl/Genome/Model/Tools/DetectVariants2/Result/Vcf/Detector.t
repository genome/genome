#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN{
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
};

use above "Genome";
use Test::More;
use Genome::Test::Factory::SoftwareResult::User;

if (Genome::Sys->arch_os ne 'x86_64') {
    plan skip_all => 'requires 64-bit machine';
}

use_ok('Genome::Model::Tools::DetectVariants2::Result::Vcf::Detector');

my $refbuild_id = 101947881;
my $ref_seq_build = Genome::Model::Build::ImportedReferenceSequence->get($refbuild_id);
ok($ref_seq_build, 'human36 reference sequence build') or die;

my $result_users = Genome::Test::Factory::SoftwareResult::User->setup_user_hash(
    reference_sequence_build => $ref_seq_build,
);

#TODO this could really use its own very tiny dataset--we don't care about the results in this test so much as the process
my $test_dir = Genome::Config::get('test_inputs') . '/Genome-Model-Tools-DetectVariants2-Samtools/';
my $test_working_dir = File::Temp::tempdir('DetectVariants2-ResultXXXXX', CLEANUP => 1, TMPDIR => 1);
my $bam_input = $test_dir . '/alignments/102922275_merged_rmdup.bam';

my $expected_dir = $test_dir . '/expected.v4/';
ok(-d $expected_dir, "expected results directory exists");

my $version = 'r613';

my $detector_parameters = '';

my %command_params = (
    reference_build_id => $refbuild_id,
    aligned_reads_input => $bam_input,
    version => $version,
    params => $detector_parameters,
    aligned_reads_sample => 'TEST',
    output_directory => $test_working_dir . '/test',
    result_users => $result_users,
);

my $command = Genome::Model::Tools::DetectVariants2::Samtools->create(%command_params);

isa_ok($command, 'Genome::Model::Tools::DetectVariants2::Samtools', 'created samtools detector');
$command->dump_status_messages(1);
ok($command->execute, 'executed samtools command');
my $result = $command->_vcf_result;
isa_ok($result, 'Genome::Model::Tools::DetectVariants2::Result::Vcf::Detector', 'generated result');

my $snvs_output = $command->output_directory."/snvs.vcf.gz";
is(readlink($snvs_output), $result->get_vcf("snvs"), 'created symlink to snvs.vcf.gz result');

my $indels_output = $command->output_directory."/indels.vcf.gz";
is(readlink($indels_output), $result->get_vcf("indels"), 'created symlink to indels.vcf.gz result');

my $snvs_result = Genome::SoftwareResult->get($result->id);

isa_ok($snvs_result, 'Genome::Model::Tools::DetectVariants2::Result::Vcf::Detector', 'Able to get detector-vcf result');
$DB::single=1;

done_testing();
