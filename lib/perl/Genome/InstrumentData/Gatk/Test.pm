package Genome::InstrumentData::Gatk::Test;

use strict;
use warnings;

use Genome;

require Genome::Utility::Test;

class Genome::InstrumentData::Gatk::Test {
    is => 'UR::Singleton',
    has_constant => [
        test_data_dir => {
            is_constant => 1,
            calculate => q| return $ENV{GENOME_TEST_INPUTS}.'/Genome-InstrumentData-Gatk'; |,
        },
        reference_build => {
            is_constant => 1,
            calculate_from => [qw/ test_data_dir /],
            calculate => q| 
                Genome::Model::Build::ImportedReferenceSequence->__define__(
                    name => 'Test Ref Build v1',
                    data_directory => $test_data_dir.'/reference',
                ); 
            |,
        },
        bam_source => {
            is_constant => 1,
            calculate_from => [qw/ test_data_dir reference_build /],
            calculate => q| 
                Genome::InstrumentData::AlignmentResult::Merged->__define__(
                    id => 9999,
                    output_dir => $test_data_dir.'/inputs',
                    reference_build => $reference_build,
                ); 
            |,
        },
        known_sites => {
            is_constant => 1,
            calculate_from => [qw/ known_sites_indel known_sites_snv /],
            calculate => q| return [ $known_sites_indel, $known_sites_snv ]; |,
        },
        known_site => { # use when running GATK
            is_constant => 1,
            calculate_from => [qw/ known_sites_indel /],
            calculate => q| return $known_sites_indel |,
        },
        known_sites_indel => {
            is_constant => 1,
            calculate_from => [qw/ indel_result /],
            calculate => q|
                Genome::Model::Build::ImportedVariationList->__define__(
                    id => 9998,
                    indel_result => $indel_result,
                );
            |,
        },
        indel_result => {
            is_constant => 1,
            calculate_from => [qw/ test_data_dir /],
            calculate => q|
                Genome::Model::Tools::DetectVariants2::Result::Manual->__define__(
                    id => 9997,
                    output_dir => $test_data_dir.'/Mills_and_1000G_gold_standard/',
                );
            |,
        },
        known_sites_snv => {
            is_constant => 1,
            calculate_from => [qw/ snv_result /],
            calculate => q|
                Genome::Model::Build::ImportedVariationList->__define__(
                    id => 9996,
                    snv_result => $snv_result,
                );
            |,
        },
        snv_result => {
            is_constant => 1,
            calculate_from => [qw/ test_data_dir /],
            calculate => q|
                Genome::Model::Tools::DetectVariants2::Result::Manual->__define__(
                    id => 9995,
                    output_dir => $test_data_dir.'/dbsnp_138/',
                );
            |,
        },
    ],
};

1;

