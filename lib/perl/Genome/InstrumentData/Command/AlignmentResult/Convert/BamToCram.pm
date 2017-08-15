package Genome::InstrumentData::Command::AlignmentResult::Convert::BamToCram;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::AlignmentResult::Convert::BamToCram {
    is => 'Genome::InstrumentData::Command::AlignmentResult::Convert::Base',
};

sub shortcut {
    my $self = shift;

    return $self->_is_currently_cram;
}

sub execute {
    my $self = shift;

    my $guard = $self->_lock()->unlock_guard;

    return 1 if $self->_is_currently_cram;

    unless($self->_is_currently_bam) {
        $self->fatal_message("This result cannot be converted to CRAM.");
    }

    my $result = $self->alignment_result;
    my $bam_file = $result->bam_file;

    my $cram_file = $bam_file;
    $cram_file =~ s/bam$/cram/;

    $self->_run_conversion('-C', $bam_file, $cram_file);

    unless ($self->_verify_file($cram_file)) {
        $self->fatal_message('Failed to convert to CRAM file.');
    } else {
        unlink $bam_file;
        $result->filetype('cram');
    }

    return 1;
}

1;
