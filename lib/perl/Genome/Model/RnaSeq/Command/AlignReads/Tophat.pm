package Genome::Model::RnaSeq::Command::AlignReads::Tophat;

use strict;
use warnings;

use Genome;

my $DEFAULT_LSF_RESOURCE = "-R 'select[model!=Opteron250 && type==LINUX64 && mem>64000 && tmp>150000] span[hosts=1] rusage[tmp=150000, mem=64000]' -M 64000000 -n 4";

class Genome::Model::RnaSeq::Command::AlignReads::Tophat {
    is => ['Command::V2'],
    has_input_output => [
        build_id => {},
    ],
    has => [
        build => { is => 'Genome::Model::Build::RnaSeq', id_by => 'build_id', },
        _unaligned_bam_files => {
            is_optional => 1,
        },
    ],
    has_param => [
        lsf_resource => { default_value => $DEFAULT_LSF_RESOURCE },
    ],
};

sub shortcut {
    my $self = shift;

    #try to get using the lock in order to wait here in shortcut if another process is creating this alignment result
    my $alignment = $self->build->alignment_result_with_lock;
    unless($alignment) {
        $self->debug_message('No existing alignment found.');
        return;
    }

    $self->_link_build_to_alignment($alignment);
    $self->debug_message('Using existing alignment ' . $alignment->__display_name__);
    return 1;
}

sub execute {
    my $self = shift;

    my $alignment = $self->build->generate_alignment_result;
    unless($alignment) {
        $self->error_message('Failed to generate alignment.');
        die $self->error_message;
    }

    $self->_link_build_to_alignment($alignment);
    $self->debug_message('Generated alignment.');
    return 1;
}

sub _link_build_to_alignment {
    my $self = shift;
    my $alignment = shift;

    my $link = $alignment->add_user(user => $self->build, label => 'uses');
    if ($link) {
        $self->debug_message("Linked alignment " . $alignment->id . " to the build");
    }
    else {
        $self->error_message(
            "Failed to link the build to the alignment "
            . $alignment->__display_name__
            . "!"
        );
        die $self->error_message;
    }

    Genome::Sys->create_symlink($alignment->output_dir, $self->build->accumulated_alignments_directory);

    return 1;
}

1;
