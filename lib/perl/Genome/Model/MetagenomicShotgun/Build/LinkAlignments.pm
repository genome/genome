package Genome::Model::MetagenomicShotgun::Build::LinkAlignments;

use strict;
use warnings;

use Genome;

class Genome::Model::MetagenomicShotgun::Build::LinkAlignments {
    is => 'Command::V2',
    has_input => [
        input_build => {
            is => 'Genome::Model::Build::MetagenomicShotgun',
            is_many => 1,
            doc => 'The MetaShot build to work with.',
        },
    ],
    has_output => [
        build => {
            is => 'Genome::Model::Build::MetagenomicShotgun',
            calculate_from => ['input_build'],
            calculate => sub{ return $_[0]; },
        },
    ],
};

sub execute {
    my $self = shift;

    for my $sub_model_label ( Genome::Model::MetagenomicShotgun->sub_model_labels ) {
        my $link_alignments_ok = $self->_link_sub_build_alignments_to_build($sub_model_label);
        return if not $link_alignments_ok;
    }

    return 1;
}

sub _link_sub_build_alignments_to_build {
    my ($self, $sub_model_label) = @_;

    Carp::confess('No sub model name given to link alignments!') if not $sub_model_label;

    my $sub_build = $self->build->model->last_complete_build_for_sub_model($sub_model_label);
    return if not $sub_build;

    my $dir = $self->build->data_directory;
    my $sub_dir = $dir.'/'.$sub_model_label;
    my $create_ok = eval{ Genome::Sys->create_directory($sub_dir); };
    if ( not $create_ok ) {
        $self->error_message($@) if $@;
        $self->error_message("Failed to create $sub_model_label sub dir! ".$sub_dir);
        return;
    }

    for my $instrument_data ( $sub_build->instrument_data ) {
        my @alignments = $sub_build->alignment_results_for_instrument_data($instrument_data); # This should only be one.
        for my $alignment ( @alignments ) {
            my $target = $alignment->output_dir;
            my $link = $sub_dir.'/'.$instrument_data->id.'-'.$alignment->id;
            unlink $link;
            my $link_ok = eval{ Genome::Sys->create_symlink($target, $link); };
            if ( not $link_ok ) {
                $self->error_message($@) if $@;
                $self->error_message("Failed to create symlink! From $link to $target");
                return;
            }
        }
    }

    return 1;
}

1;

