package Genome::Qc::Tool::VerifyBamId;

use strict;
use warnings;
use Genome;
use List::MoreUtils qw(uniq);

class Genome::Qc::Tool::VerifyBamId {
    is => 'Genome::Qc::Tool::WithVariationListVcf',
};


sub cmd_line {
    my $self = shift;

    my $cmd = $self->gmt_class->create($self->gmt_params);
    return $cmd->_get_cmd_list;
}

sub metrics {
    return (
        number_snps => '#SNPS',
        freemix => 'FREEMIX',
        chipmix => 'CHIPMIX',
    );
}

sub get_metrics {
    my $self = shift;

    my $file = $self->qc_metrics_file . '.selfSM';
    my $reader = Genome::Utility::IO::SeparatedValueReader->create(
        input => $file,
        separator => '\t',
        is_regex => 1,
    );

    my %metrics = $self->metrics;
    while ( my $line = $reader->next ) {
        if ($line->{'#SEQ_ID'} eq $self->sample_name) {
            my %desired_metric_results;
            while (my ($metric_name, $metric_key) = each %metrics) {
                $desired_metric_results{$metric_name} = $line->{$metric_key};
            }
            return %desired_metric_results;
        }
    }
    return ();
}

sub gmt_class {
    return 'Genome::Model::Tools::VerifyBamId';
}

sub qc_metrics_file_accessor {
    return 'out_prefix';
}

1;
