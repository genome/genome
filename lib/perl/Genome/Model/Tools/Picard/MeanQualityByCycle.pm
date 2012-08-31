package Genome::Model::Tools::Picard::MeanQualityByCycle;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Picard::MeanQualityByCycle {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        input_file   => {
            is  => 'String',
            doc => 'The SAM/BAM files to run on.  File type is determined by suffix.',
        },
        output_file  => {
            is  => 'String',
            doc => 'File to write quality score metrics to.',
        },
        chart_output  => {
            is  => 'String',
            doc => 'File to write chart to.',
        },
        reference_sequence => {
            is_optional => 1,
        },
        stop_after => {
            is  => 'Integer',
            doc => 'Stop after processing N reads, mainly for debugging. Default value: 0.',
            is_optional   => 1,
        },
        aligned_reads_only => {
            is => 'Boolean',
            default_value => 0,
            is_optional => 1,
        },
        pf_reads_only => {
            is => 'Boolean',
            default_value => 0,
            is_optional => 1,
        },
        assume_sorted => {
            is => 'Boolean',
            default_value => 1,
            is_optional => 1,
        },
    ],
};

sub help_brief {
    'Reads a SAM or BAM file and writes a file containing metrics about the statistical distribution of insert size (excluding duplicates) and generates a histogram plot.';
}

sub help_detail {
    return <<EOS
    For Picard documentation of this command see:
    http://picard.sourceforge.net/command-line-overview.shtml#CollectInsertSizeMetrics
EOS
}

sub execute {
    my $self = shift;

    my $cmd = $self->picard_path .'/MeanQualityByCycle.jar net.sf.picard.analysis.MeanQualityByCycle';
    $cmd   .= ' OUTPUT='. $self->output_file  .' INPUT='. $self->input_file .' CHART_OUTPUT='. $self->chart_output;
    if (defined($self->stop_after)) {
        $cmd .= ' STOP_AFTER='. $self->stop_after;
    }
    if (defined($self->reference_sequence)) {
        $cmd .= ' REFERENCE_SEQUENCE='. $self->reference_sequence;
    }
    if (defined($self->aligned_reads_only)) {
        if ($self->aligned_reads_only) {
            $cmd .= ' ALIGNED_READS_ONLY=true';
        } else {
            $cmd .= ' ALIGNED_READS_ONLY=false';
        }
    }
    if (defined($self->pf_reads_only)) {
        if ($self->pf_reads_only) {
            $cmd .= ' PF_READS_ONLY=true';
        } else {
            $cmd .= ' PF_READS_ONLY=false';
        }
    }
    if (defined($self->assume_sorted)) {
        if ($self->assume_sorted) {
            $cmd .= ' ASSUME_SORTED=true';
        } else {
            $cmd .= ' ASSUME_SORTED=false';
        }
    }
    $self->run_java_vm(
        cmd          => $cmd,
        input_files  => [$self->input_file],
        output_files => [$self->output_file, $self->chart_output],
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
