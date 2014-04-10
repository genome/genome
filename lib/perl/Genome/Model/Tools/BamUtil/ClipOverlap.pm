package Genome::Model::Tools::BamUtil::ClipOverlap;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::BamUtil::ClipOverlap {
    doc => "Run BamUtil with the 'ClipOverlap' tool",
    is => 'Genome::Model::Tools::BamUtil',
    has_input => [
        input_file => {
            is => 'Text',
            doc => 'The path to the original bam/sam you would like to be clipped',
        },
        output_file => {
            is => 'Text',
            doc => "The path to the clipped bam/sam",
        },
        index_bam => {
            is => 'Boolean',
            default => 1,
            doc => 'Index the bam after alignment.'
        },
        file_format => {
            is => 'Text',
            default => 'bam',
            valid_values => ['sam', 'bam'],
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

    if ($self->index_bam && $self->file_format eq 'bam') {
        my $clipped_bam = $self->output_file;
        die $self->error_message("Couldn't find clipped bam at $clipped_bam!") unless -f $clipped_bam;

        my $rv = Genome::Model::Tools::Sam::IndexBam->execute(bam_file => $clipped_bam);
        die $self->error_message("Failed to run gmt sam index-bam on $clipped_bam") unless $rv->result == 1;
    }

    return 1;
}

sub clipoverlap_creator_command {
    my $self = shift;
    my $bamutil_command = $self->base_command;
    $bamutil_command .= " clipoverlap";
    $bamutil_command .= " --in " . $self->input_file;
    $bamutil_command .= " --out " . $self->output_file;
    return $bamutil_command;
}

sub _check_inputs {
    my $self = shift;

    Genome::Sys->validate_file_for_reading($self->input_file);
    Genome::Sys->validate_file_for_writing($self->output_file);

    return 1;
}

1;
