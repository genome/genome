package Genome::Model::MetagenomicShotgun::Build::AlignToMetaNr;

use strict;
use warnings;

use Genome;

class Genome::Model::MetagenomicShotgun::AlignToMetaNr {
    is => 'Genome::Model::MetagenomicShotgun::Base',
    has_ouput => [
        aligned => { is_many => 1, },
    ],
};

sub execute {
    my $self = shift;

    my $build = $self->build;
    my $metagenomic_protein_model = $build->model->metagenomic_protein_model;
    my $mg_protein_build = $self->_start_build($metagenomic_protein_model, $self->instrument_data);
    my $mg_nr_build_ok = $self->_wait_for_build($mg_protein_build);
    return if not $mg_nr_build_ok;

    my $link_alignments = $self->_link_sub_build_alignments_to_build(build => $build, sub_build => $mg_protein_build, sub_model_name => 'metagenomic_protein');
    return if not $link_alignments;

    my @mg_protein_aligned = $self->_extract_data->($mg_protein_build, "aligned");
    $self->aligned(@mg_protein_aligned);

    return 1;
}

1;

