package Genome::Annotation::Readcount;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Readcount {
    is => 'Genome::Annotation::Detail::Command',
    has_input => [
        aligned_bam_result => {
            is => 'Genome::InstrumentData::AlignedBamResult',
            doc => 'The bam result used to calculate read counts',
        },
        vcf_result => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Vcf',
            doc => "The vcf result that will provide positions where read count information will be generated",
        },
        variant_type => {
            is => 'Text',
            doc => "The type of variant the vcf_result represents",
            valid_values => ['snvs', 'indels'],
        },
        use_version => {
            is  => 'Version',
            doc => "bam-readcount version to be used.",
        },
        minimum_mapping_quality => {
            is => 'Integer',
            default => 0,
            doc => "filter reads with mapping quality less than this. This is the -q parameter.",
        },
        minimum_base_quality => {
            is => 'Integer',
            default => 0,
            doc => "don't include reads where the base quality is less than this. This is the -b parameter. This is only available in versions 0.3 and later.",
        },
        max_count => {
            is  => 'Integer',
            default => 10_000_000,
            doc => "max depth to avoid excessive memory. This is the -d parameter in version 0.5.",
        },
        per_library => {
            is  => 'Bool',
            default => 0,
            doc => "report results per library. This is the -p parameter in version 0.5.",
        },
        insertion_centric => {
            is  => 'Bool',
            default => 0,
            doc => "do not include reads containing insertions after the current position in per-base counts. This is the -i parameter in version 0.5.",
        },
    ],
    has_optional_output => [
        software_result => {
            is => 'Genome::Annotation::Readcount::Result',
            doc => 'The software result created during command execution',
        },
    ],
};

sub execute {
    my $self = shift;
    die "You must supply a version greater than or equal to 0.5" unless $self->use_version >= 0.5;

    $self->software_result(Genome::Annotation::Readcount::Result->get_or_create($self->input_hash));
    return 1;
}
