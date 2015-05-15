package Genome::Qc::Tool;

use strict;
use warnings;
use Genome;

class Genome::Qc::Tool {
    is_abstract => 1,
    has => [
        alignment_result => {
            is => 'Genome::InstrumentData::AlignedBamResult',
        },
        gmt_params => {
            is => 'Hash',
        }
    ],
    has_calculated => [
        bam_file => {
            is => 'Path',
            calculate_from => [qw(alignment_result)],
            calculate => q{return $alignment_result->get_bam_file},
        },
    ],
};

sub cmd_line {
    my $self = shift;
    die $self->error_message("Abstract method run must be overriden by subclass");
}

sub output_file {
    my $self = shift;
    my $output_file_accessor = $self->output_file_accessor;
    my %params = %{$self->gmt_params};
    return $params{$output_file_accessor};
}

sub supports_streaming {
    my $self = shift;
    die $self->error_message("Abstract method supports_streaming must be overridden by subclass");
}

sub get_metrics {
    my $self = shift;
    die $self->error_message("Abstract method get_metrics must be overridden by subclass");
}

# Overwrite this in subclass to return the gmt tool parameter name for the output file
sub output_file_accessor {
    return undef;
}

sub reference_build {
    my $self = shift;
    return $self->alignment_result->reference_build;
}

sub reference_sequence {
    my $self = shift;
    return $self->reference_build->full_consensus_path('fa');
}

1;
