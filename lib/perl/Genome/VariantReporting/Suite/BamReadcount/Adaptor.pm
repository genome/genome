package Genome::VariantReporting::Suite::BamReadcount::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Suite::BamReadcount::Adaptor {
    is => 'Genome::VariantReporting::Framework::Component::Adaptor',

    has_planned_output => [
        aligned_bam_result_id => {
            is_many => 1,
            is => 'Text',
            is_translated => 1,
            doc => 'The bam result used to calculate read counts',
        },
        version => {
            is  => 'Version',
            doc => "bam-readcount version to be used.",
        },
        minimum_mapping_quality => {
            is => 'Integer',
            doc => "filter reads with mapping quality less than this. This is the -q parameter.",
        },
        minimum_base_quality => {
            is => 'Integer',
            doc => "don't include reads where the base quality is less than this. This is the -b parameter. This is only available in versions 0.3 and later.",
        },
        max_count => {
            is  => 'Integer',
            doc => "max depth to avoid excessive memory. This is the -d parameter in version 0.5.",
        },
        per_library => {
            is  => 'Boolean',
            doc => "report results per library. This is the -p parameter in version 0.5.",
        },
        insertion_centric => {
            is  => 'Boolean',
            doc => "do not include reads containing insertions after the current position in per-base counts. This is the -i parameter in version 0.5.",
        },
    ],
    doc => 'Add bam readcount information to a vcf',
};

sub name {
    "bam-readcount";
}

1;
