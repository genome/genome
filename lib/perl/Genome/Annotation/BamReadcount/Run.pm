package Genome::Annotation::BamReadcount::Run;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::BamReadcount::Run {
    is => 'Genome::Annotation::CommandBase',
    has_input => [
        aligned_bam_result => {
            is => 'Genome::InstrumentData::AlignedBamResult',
            doc => 'The bam result used to calculate read counts',
        },
        version => {
            is  => 'Version',
            example_values => ['0.5'],
            doc => "bam-readcount version to be used.",
        },
        minimum_mapping_quality => {
            is => 'Integer',
            example_values => [0],
            doc => "filter reads with mapping quality less than this. This is the -q parameter.",
        },
        minimum_base_quality => {
            is => 'Integer',
            example_values => [0],
            doc => "don't include reads where the base quality is less than this. This is the -b parameter. This is only available in versions 0.3 and later.",
        },
        max_count => {
            is  => 'Integer',
            example_values => [10_000_000],
            doc => "max depth to avoid excessive memory. This is the -d parameter in version 0.5.",
        },
        per_library => {
            is  => 'Bool',
            example_values => [0],
            doc => "report results per library. This is the -p parameter in version 0.5.",
        },
        insertion_centric => {
            is  => 'Bool',
            example_values => [0],
            doc => "do not include reads containing insertions after the current position in per-base counts. This is the -i parameter in version 0.5.",
        },
    ],
};

sub execute {
    my $self = shift;
    die "You must supply a version greater than or equal to 0.5" unless $self->version >= 0.5;

    $self->output_result(Genome::Annotation::BamReadcount::RunResult->get_or_create($self->input_hash));
    return 1;
}
