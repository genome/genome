package Genome::Qc::Tool::Samtools::Flagstat;

use strict;
use warnings;
use Genome;

class Genome::Qc::Tool::Samtools::Flagstat {
    is => 'Genome::Qc::Tool',
};

sub supports_streaming {
    return 1;
}

sub cmd_line {
    my $self = shift;
    my @cmd = qw(gmt sam flagstat);
    while (my ($param_name, $param_value) = each %{$self->gmt_params}) {
        push @cmd, "--$param_name=$param_value";
    }
    return @cmd;
}

sub get_metrics {
    my $self = shift;

    my $file = $self->qc_metrics_file;
    return %{Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($file)};
}

sub qc_metrics_file_accessor {
    return 'output-file';
}

1;
