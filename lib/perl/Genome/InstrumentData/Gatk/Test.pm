package Genome::InstrumentData::Gatk::Test;

use strict;
use warnings;

use Genome;

require Genome::Utility::Test;

class Genome::InstrumentData::Gatk::Test {
    is => 'UR::Singleton',
    has_constant => [
        test_data_dir => {
            calculate => q| return $ENV{GENOME_TEST_INPUTS}.'/Genome-InstrumentData-Gatk'; |,
        },
        reference_build => {
            calculate_from => [qw/ test_data_dir /],
            calculate => q| 
                Genome::Model::Build::ImportedReferenceSequence->__define__(
                    name => 'Test Ref Build v1',
                    data_directory => $test_data_dir.'/reference',
                ); 
            |,
        },
        bam_source => {
            calculate_from => [qw/ test_data_dir /],
            calculate => q| 
                Genome::InstrumentData::AlignmentResult::Merged->__define__(
                    id => 9999,
                    output_dir => $test_data_dir.'/inputs',
                ); 
            |,
        },
        known_indel => {
            calculate_from => [qw/ indel_result /],
            calculate => q|
                Genome::Model::Build::ImportedVariationList->__define__(
                    id => 9998,
                    indel_result => $indel_result,
                );
            |,
        },
        indel_result => {
            calculate_from => [qw/ test_data_dir /],
            calculate => q|
                Genome::Model::Tools::DetectVariants2::Result::Manual->__define__(
                    id => 9997,
                    output_dir => $test_data_dir.'/Mills_and_1000G_gold_standard/',
                );
            |,
        },
    ],
};

1;

