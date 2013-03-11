package Genome::Model::MetagenomicShotgun::Build::AlignTo;

use strict;
use warnings;

use Genome;

class Genome::Model::MetagenomicShotgun::AlignTo{
    is => 'Genome::Model::MetagenomicShotgun::Base',
    has_param => [
        sub_model_label => { 
            is => 'Text',
            valid_values => [ Genome::Model::MetagenomicShotgun->sub_model_labels ],
        },
    ],
    has_opitonal_ouput => [
        aligned => { is_many => 1, },
        unaligned => { is_many => 1, },
    ],
};

sub execute {
    my $self = shift;

    my $build = $self->build;
    my $sub_model_label = $self->sub_model_label;
    my $viral_nucleotide_model = $build->model->$sub_model_label;
    my @instrument_data = $self->instrument_data;
    my $viral_nucleotide_build = $self->_start_build($viral_nucleotide_model, @instrument_data);
    my $viral_nt_build_ok = $self->_wait_for_build($viral_nucleotide_build);
    return if not $viral_nt_build_ok;

    my $link_alignments = $self->_link_sub_build_alignments_to_build(
        build => $build,
        sub_build => $viral_nucleotide_build,
        sub_model_name => $self->sub_model_label,
    );
    return if not $link_alignments;

    return 1;
}

1;

