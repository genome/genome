package Genome::Model::MetagenomicShotgun::Build::AlignToMetaNt;

use strict;
use warnings;

use Genome;

class Genome::Model::MetagenomicShotgun::AlignToMetaNt {
    is => 'Genome::Model::MetagenomicShotgun::Base',
    has_ouput => [
        instrument_data => { is_many => 1, },
    ],
    has_ouput => [
        aligned => { is_many => 1, },
        unaligned => { is_many => 1, },
    ],
};

sub execute {
    my $self = shift;

    my $build = $self->build;
    my $metagenomic_nucleotide_model = $build->model->metagenomic_nucleotide_model;
    my $mg_nucleotide_build = $self->_start_build($metagenomic_nucleotide_model, $self->instrument_data);
    my $mg_nt_build_ok = $self->_wait_for_build($mg_nucleotide_build);
    return if not $mg_nt_build_ok;

    my $link_alignments = $self->_link_sub_build_alignments_to_build(build => $build, sub_build => $mg_nucleotide_build, sub_model_name => 'metagenomic_nucleotide');
    return if not $link_alignments;

    my @mg_nucleotide_aligned = $self->_extract_data($mg_nucleotide_build, "aligned");
    $self->aligned(@mg_nucleotide_aligned);

    my @mg_nucleotide_unaligned = $self->_extract_data($mg_nucleotide_build, "unaligned");
    $self->unaligned(@mg_nucleotide_unaligned);

    return 1;
}

1;

