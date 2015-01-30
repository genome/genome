package Genome::Model::Tools::Vcf::CreateCrossSampleVcf::Process;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Vcf::CreateCrossSampleVcf::Process {
    is => 'Genome::Process',
    has_input => [
        builds => {
            is => 'Genome::Model::Build',
            is_many => 1,
        },
        roi_list => {
            is => 'Genome::FeatureList',
            is_optional => 1,
        },
        wingspan => {
            is => 'Text',
        },
        joinx_version => {
            is => 'Text',
        },
    ],
};

1;
