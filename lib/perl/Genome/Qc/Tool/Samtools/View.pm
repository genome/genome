package Genome::Qc::Tool::Samtools::View;

use strict;
use warnings;
use Genome;

class Genome::Qc::Tool::Samtools::View {
    is => 'Genome::Qc::Tool',
};

sub supports_streaming {
    return 1;
}

sub cmd_line {
    my $self = shift;
    return (qw(samtools view), $self->gmt_params->{input_file});
}

sub get_metrics {
    return {};
}

1;

