package Genome::Qc::Tool;

use strict;
use warnings;
use Genome;

class Genome::Qc::Tool {
    has => [
        alignment_result => {
            is => 'Genome::InstrumentData::AlignmentResult',
        },
        gmt_params => {
            is => 'Hash',
        }
    ],
    has_calculated => [
        bam_file => {
            is => 'Path',
            calculate_from => [qw(alignment_result)],
            calculate => q{return $alignment_result->bam_path},
        },
    ],
};

sub cmd_line {
    my $self = shift;
    die $self->error_message("Abstract method run must be overriden by subclass");
}

sub supports_streaming {
    my $self = shift;
    die $self->error_message("Abstract method supports_streaming must be overridden by subclass");
}

sub get_metrics {
    my $self = shift;
    die $self->error_message("Abstract method get_metrics must be overridden by subclass");
}

1;

