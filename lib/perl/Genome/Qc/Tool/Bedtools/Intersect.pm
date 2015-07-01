package Genome::Qc::Tool::Bedtools::Intersect;

use strict;
use warnings;
use Genome;

class Genome::Qc::Tool::Bedtools::Intersect {
    is => 'Genome::Qc::Tool::WithVariationListVcf',
};

sub cmd_line {
    my $self = shift;

    my @cmd = qw(gmt bed-tools intersect);
    while (my ($input, $value) = each %{$self->gmt_params}) {
        if ($input eq 'header') {
            push @cmd, "--$input";
        }
        else {
            push @cmd, "--$input=$value";
        }
    }
    return @cmd;
}

sub get_metrics {
    return ();
}

sub target_region_bed_file {
    my $self = shift;
    return $self->alignment_result->instrument_data->target_region_set->processed_bed_file;
}

1;
