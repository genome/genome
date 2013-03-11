package Genome::Model::MetagenomicShotgun::Build::ScreenContamination;

use strict;
use warnings;

use Genome;

class Genome::Model::MetagenomicShotgun::ScreenContamination{
    is => 'Genome::Model::MetagenomicShotgun::Base',
    has_ouput => [
        unaligned => {
            is_many => 1,
        },
    ],
};

sub execute {
    my $self = shift;

    my $build = $self->build;
    my $contamination_screen_model = $build->model->contamination_screen_model;
    my @original_instdata = $build->instrument_data;
    my $cs_build = $self->_start_build($contamination_screen_model, @original_instdata);
    my $cs_build_ok = $self->_wait_for_build($cs_build);
    return if not $cs_build_ok;
    my $link_alignments = $self->_link_sub_build_alignments_to_build(build => $build, sub_build => $cs_build, sub_model_name => 'contamination_screen');
    return if not $link_alignments;

    my @cs_unaligned;
    if ($build->model->filter_contaminant_fragments){
        @cs_unaligned = $self->_extract_data($cs_build, "unaligned paired");
    }
    else{
        @cs_unaligned = $self->_extract_data($cs_build, "unaligned");
    }

    $self->unaligned(@cs_unaligned);

    return 1;
}

1;

