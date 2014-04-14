package Genome::Annotation::AddReadcounts;

use strict;
use warnings FATAL => 'all';

class Genome::Annotation::AddReadcounts {
    is => 'Genome::Annotation::Detail::Command',
    has_input => [
        input_vcf_result => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Vcf',
        },
        variant_type => {
            is => 'Text',
            valid_values => ['snvs', 'indels'],
        },
        readcount_results => {
            is => 'Genome::Annotation::Readcount::Result',
            is_many => 1,
        },
    ],
    has_output => [
        software_result => {
            is => 'Genome::Annotation::AddReadcounts::Result',
        },
    ],
};

sub execute {
    my $self = shift;

    $self->software_result(Genome::Annotation::AddReadcounts::Result->get_or_create($self->input_hash));
    return 1;
}

