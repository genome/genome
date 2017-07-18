package Genome::InstrumentData::AlignmentResult::Command::CramToBam;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::AlignmentResult::Command::CramToBam {
    is => 'Genome::InstrumentData::AlignmentResult::Command::ConvertBase',
};

sub shortcut {
    my $self = shift;

    return $self->_is_currently_bam;
}

sub execute {
    my $self = shift;

    my $guard = $self->_lock()->unlock_guard;

    return 1 if $self->_is_currently_bam;

    unless($self->_is_currently_cram) {
        $self->fatal_message("This result cannot be converted to BAM.");
    }

    my $result = $self->alignment_result;
    my $cram_file = $result->bam_file;

    my $bam_file = $cram_file;
    $bam_file =~ s/cram$/bam/;

    $self->_run_conversion('-b', $cram_file, $bam_file);

    unless ($self->_verify_file($bam_file)) {
        $self->fatal_message('Failed to convert to BAM file.');
    } else {
        unlink $cram_file;
        $result->filetype('bam');
    }

    return 1;
}

1;
