#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use File::Temp;
use Test::More tests => 20;
use Data::Dumper;
use File::Compare;
use Genome::Test::Factory::SoftwareResult::User;

$ENV{UR_DBI_NO_COMMIT} = 1;
$ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;

BEGIN {
    use_ok( 'Genome::Model::Tools::DetectVariants2::Filter::VarscanHighConfidence')
};

my $refbuild_id = 101947881;
my $test_data_directory = Genome::Config::get('test_inputs') . "/Genome-Model-Tools-DetectVariants2-Filter-VarscanHighConfidence";

my $result_users = Genome::Test::Factory::SoftwareResult::User->setup_user_hash(
    reference_sequence_build_id => $refbuild_id,
);

# Updated to .v2 for correcting an error with newlines
my $expected_directory = $test_data_directory . "/expected.v3";
my $detector_directory = $test_data_directory . "/varscan-somatic-2.2.4-";
my $tumor_bam_file  = $test_data_directory. "/alignments/102922275_merged_rmdup.bam"; #'/flank_tumor_sorted.bam';
my $normal_bam_file  = $tumor_bam_file;
my $test_output_base = File::Temp::tempdir('Genome-Model-Tools-DetectVariants2-Filter-VarscanHighConfidence-XXXXX', CLEANUP => 1, TMPDIR => 1);
my $test_output_dir = $test_output_base . '/filter';

my $hq_output_bed = "$test_output_dir/indels.hq.bed";
my $lq_output_bed = "$test_output_dir/indels.lq.bed";

my @expected_output_files = qw| snvs.Germline     
                                snvs.Germline.lc  
                                snvs.LOH.hc  
                                snvs.Somatic     
                                snvs.Somatic.lc  
                                snvs.hq.bed     
                                snvs.hq.v2.bed  
                                snvs.lq.bed     
                                snvs.lq.v2.bed
                                snvs.Germline.hc  
                                snvs.LOH          
                                snvs.LOH.lc  
                                snvs.Somatic.hc  
                                snvs.hq          
                                snvs.hq.v1.bed  
                                snvs.lq         
                                snvs.lq.v1.bed |;

my $detector_result = Genome::Model::Tools::DetectVariants2::Result->__define__(
    output_dir => $detector_directory,
    detector_name => 'test',
    detector_params => '',
    detector_version => 'awesome',
    aligned_reads => $tumor_bam_file,
    control_aligned_reads => $normal_bam_file,
    reference_build_id => $refbuild_id,
);
$detector_result->lookup_hash($detector_result->calculate_lookup_hash);

my $varscan_high_confidence = Genome::Model::Tools::DetectVariants2::Filter::VarscanHighConfidence->create(
    previous_result_id => $detector_result->id,
    output_directory => $test_output_dir,
    result_users => $result_users,
);

ok($varscan_high_confidence, "created VarscanHighConfidence object");
ok($varscan_high_confidence->execute(), "executed VarscanHighConfidence");

my @list = glob($test_output_dir."/*");

for my $output_file (@expected_output_files){
    my $expected_file = $expected_directory."/".$output_file;
    my $actual_file = $test_output_dir."/".$output_file;
    is(compare($actual_file, $expected_file), 0, "$actual_file output matched expected output");
}
