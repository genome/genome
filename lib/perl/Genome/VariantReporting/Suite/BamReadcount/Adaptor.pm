package Genome::VariantReporting::Suite::BamReadcount::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;
use version 0.77;
use Params::Validate qw(validate_pos);

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
            doc => "don't include reads where the base quality is less than this. This is the -b parameter.",
        },
        max_count => {
            is  => 'Integer',
            doc => "max depth to avoid excessive memory. This is the -d parameter.",
        },
        per_library => {
            is  => 'Boolean',
            doc => "report results per library. This is the -p parameter.",
        },
        insertion_centric => {
            is  => 'Boolean',
            doc => "do not include reads containing insertions after the current position in per-base counts. This is the -i parameter.",
        },
    ],
    doc => 'Add bam readcount information to a vcf',
};

my $MIN_VERSION = '0.7';

sub __planned_output_errors__ {
    my ($self, $params) = validate_pos(@_, 1, 1);

    my @errors = $self->SUPER::__planned_output_errors__($params);

    if (version->parse("v".$params->{version}) < version->parse("v$MIN_VERSION")) {
        push @errors, UR::Object::Tag->create(
            type => 'error',
            properties => ['version'],
            desc => sprintf("The BamReadcount expert requires version (%s) ".
                    "or higher, but you supplied version (%s).",
                    $MIN_VERSION, $params->{version} || 'undef'),
        );
    }
    return @errors;
}

sub name {
    "bam-readcount";
}

1;
