package Genome::Qc::Tool::Picard;

use strict;
use warnings;
use Genome;

use Hash::Flatten;

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
    my $metric_results = $gmt_class->_parse_metrics_file_into_hashref($file, undef, undef, $self->gmt_params->{metric_accumulation_level});
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
    my ($self, $picard_metrics_hash) = @_;
    my %flattened_metrics_hash;

    my $flattener = Hash::Flatten->new({
        HashDelimiter => "\t",
        OnRefScalar => 'die',
        OnRefRef => 'die',
        OnRefGlob => 'die',
        ArrayDelimiter => "ERROR!",
    });

    my $flat_hash = $flattener->flatten($picard_metrics_hash);
    for my $key (keys %$flat_hash) {
        if ($key =~ /\t{2,}/) {
            my $value = delete $flat_hash->{$key};
            $key =~ s/\t{2,}/\t/g;
            $flat_hash->{$key} = $value;
        }
    }

    return %$flat_hash;
}

1;
