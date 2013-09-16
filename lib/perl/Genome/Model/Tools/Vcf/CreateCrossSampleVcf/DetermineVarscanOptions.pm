package Genome::Model::Tools::Vcf::CreateCrossSampleVcf::DetermineVarscanOptions;

use strict;
use warnings;

use above 'Genome';
use File::Spec;

class Genome::Model::Tools::Vcf::CreateCrossSampleVcf::DetermineVarscanOptions {
    is => 'Command::V2',
    has_input => [
        output_directory => {
            is => 'Text',
        },
    ],
    has_calculated_output => [
        output_file => {
            is => 'Path',
            calculate =>  q{ File::Spec->join($output_directory, "varscan_consensus.vcf") },
            calculate_from => ['output_directory'],
        },
        output_vcf => {
            is => 'Boolean',
            calculate =>  q{ 1; },
        },
    ],
};
