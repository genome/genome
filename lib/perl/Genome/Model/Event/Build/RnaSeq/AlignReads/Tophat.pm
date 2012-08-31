package Genome::Model::Event::Build::RnaSeq::AlignReads::Tophat;

use strict;
use warnings;

use version;
use Genome;

class Genome::Model::Event::Build::RnaSeq::AlignReads::Tophat {
    is => ['Genome::Model::Event::Build::RnaSeq::AlignReads'],
    has => [
        _unaligned_bam_files => {
            is_optional => 1,
        },
    ],
};

sub bsub_rusage {
    return "-R 'select[model!=Opteron250 && type==LINUX64 && mem>16000 && tmp>150000] span[hosts=1] rusage[tmp=150000, mem=16000]' -M 16000000 -n 4";
}

sub shortcut {
    my $self = shift;

    #try to get using the lock in order to wait here in shortcut if another process is creating this alignment result
    my $alignment = $self->build->alignment_result_with_lock;
    unless($alignment) {
        $self->status_message('No existing alignment found.');
        return;
    }

    $self->_link_build_to_alignment($alignment);
    $self->status_message('Using existing alignment ' . $alignment->__display_name__);
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
    $self->status_message('Generated alignment.');
    return 1;
}

sub _link_build_to_alignment {
    my $self = shift;
    my $alignment = shift;

    my $link = $alignment->add_user(user => $self->build, label => 'uses');
    if ($link) {
        $self->status_message("Linked alignment " . $alignment->id . " to the build");
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
