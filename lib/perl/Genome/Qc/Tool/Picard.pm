package Genome::Qc::Tool::Picard;

use strict;
use warnings;
use Genome;

class Genome::Qc::Tool::Picard {
    is => 'Genome::Qc::Tool',
    is_abstract => 1,
};

sub cmd_line {
    my $self = shift;
    my $cmd = $self->gmt_class->create($self->gmt_params);
    return $cmd->build_cmdline_list;
}

sub get_metrics {
    my $self = shift;

    my $file = $self->qc_metrics_file;
    my $gmt_class = $self->gmt_class;
    my $metric_results = $gmt_class->parse_file_into_metrics_hashref($file, undef, $self->gmt_params->{metric_accumulation_level});
    my $metric_header_as_key = $gmt_class->can('_metric_header_as_key') ? $gmt_class->_metric_header_as_key : undef;
    return $self->_flatten_metrics_hash($metric_results, metric_header_as_key => $metric_header_as_key);
}

sub metrics {
    my $self = shift;
    die $self->error_message("Abstract method metrics must be overridden by subclass");
}

sub gmt_class {
    my $self = shift;
    die $self->error_message("Abstract method gmt_class must be overridden by subclass");
}

1;
