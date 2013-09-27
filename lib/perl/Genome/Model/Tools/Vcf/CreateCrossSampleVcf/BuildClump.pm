package Genome::Model::Tools::Vcf::CreateCrossSampleVcf::BuildClump;

use strict;
use warnings;

use UR;


class Genome::Model::Tools::Vcf::CreateCrossSampleVcf::BuildClump {
    is => 'UR::Value::JSON',
    id_by =>  [qw(
        backfilled_vcf
        bam_file
        pileup_output_file
        sample
        vcf_file
    )],
};


1;
