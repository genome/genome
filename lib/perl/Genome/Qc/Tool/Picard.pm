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
    my $metric_results = $gmt_class->parse_file_into_metrics_hashref($file);
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

sub _flatten_metrics_hash {
    my ($self, $picard_metrics_hash, %params) = @_;
    my %flattened_metrics_hash;
    my $metric_header_as_key = $params{metric_header_as_key};

    while (my ($metric_category, $nested_metrics_hash) = each %{$picard_metrics_hash}) {
        while (my ($metric_name, $metric_value) = each %{$nested_metrics_hash}) {
            next if $metric_header_as_key eq $metric_name;
            my $key;
            if ($metric_header_as_key) {
                $key = join('-', $metric_category, $metric_name);
                $key =~ s/^$metric_header_as_key-//;
            } else {
                $key = $metric_name;
            }
            $flattened_metrics_hash{$key} = $metric_value;
        }
    }
    return %flattened_metrics_hash;
}

1;
