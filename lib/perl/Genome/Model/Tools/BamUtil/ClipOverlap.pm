package Genome::Model::Tools::BamUtil::ClipOverlap;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::BamUtil::ClipOverlap {
    doc => "Run BamUtil with the 'ClipOverlap' tool",
    is => 'Genome::Model::Tools::BamUtil',
    has_input => [
        known => {
            is => 'Text',
            doc => 'The file of known indels',
            is_optional => 1,
            is_many => 1,
        },
        input_bam => {
            is => 'Text',
            doc => 'The path to the original bam you would like to be clipped',
        },
        output_bam => {
            is => 'Text',
            doc => "The path to the clipped bam",
        },
    ],
};

sub help_brief {
    "Run BamUtil with the 'ClipOverlap' tool"
}

sub help_synopsis {
    return <<EOS
    gmt bam-util clip-overlap --in in.bam --out out.bam
EOS
}

sub execute {
    my $self = shift;

    unless ($self->_check_inputs) {
        return;
    }
    my $command = $self->clipoverlap_creator_command;

    unless (Genome::Sys->shellcmd(cmd => $command)) {
        die $self->error_message("Failed to execute $command");
    }

    return 1;
}

sub clipoverlap_creator_command {
    my $self = shift;
    my $bamutil_command = $self->base_command;
    $bamutil_command .= " clipoverlap";
    $bamutil_command .= " --in " . $self->input_bam;
    $bamutil_command .= " --out " . $self->output_bam;
    return $bamutil_command;
}

sub _check_inputs {
    my $self = shift;

    Genome::Sys->validate_file_for_reading($self->input_bam);
    Genome::Sys->validate_file_for_writing($self->output_bam);

    return 1;
}

1;
