package Genome::Model::Tools::Picard::CollectInsertSizeMetrics;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Picard::CollectInsertSizeMetrics {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        input_file   => {
            is  => 'String',
            doc => 'The SAM/BAM files to run on.  File type is determined by suffix.',
        },
        output_file  => {
            is  => 'String',
            doc => 'File to write insert size metrics to.',
        },
        histogram_file  => {
            is  => 'String',
            doc => 'File to write insert size histogram chart to.',
        },
        tail_limit => {
            is  => 'Integer',
            doc => 'When calculating mean and stdev stop when the bins in the tail of the distribution contain fewer than mode/TAIL_LIMIT items. This also limits how much data goes into each data category of the histogram. Default value: 10000.',
            is_optional => 1,
        },
        histogram_width => {
            is  => 'Integer',
            doc => 'Explicitly sets the histogram width, overriding the TAIL_LIMIT option. Also, when calculating mean and stdev, only bins <= HISTOGRAM_WIDTH will be included.Default value: null.',
            is_optional => 1,
        },
        minimum_pct  => {
            is  => 'Integer',
            doc => 'When generating the histogram, discard any data categories (out of FR, TANDEM, RF) that have fewer than this percentage of overall reads. (Range: 0 to 1) Default value: 0.01.',
            is_optional   => 1,
        },
        stop_after => {
            is  => 'Integer',
            doc => 'Stop after processing N reads, mainly for debugging. Default value: 0.',
            is_optional   => 1,
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

    my $cmd = $self->picard_path .'/CollectInsertSizeMetrics.jar net.sf.picard.analysis.CollectInsertSizeMetrics';
    $cmd   .= ' OUTPUT='. $self->output_file  .' INPUT='. $self->input_file .' HISTOGRAM_FILE='. $self->histogram_file;

    if (defined($self->tail_limit)) {
        $cmd .= ' TAIL_LIMIT=' . $self->tail_limit;
    }
    if (defined($self->histogram_width)) {
        $cmd .= ' HISTOGRAM_WIDTH='. $self->histogram_width;
    }
    if (defined($self->minimum_pct)) {
        $cmd .= ' MINIMUM_PCT='. $self->minimum_pct;
    }
    if (defined($self->stop_after)) {
        $cmd .= ' STOP_AFTER='. $self->stop_after;
    }
    $self->run_java_vm(
        cmd          => $cmd,
        input_files  => [$self->input_file],
        output_files => [$self->output_file, $self->histogram_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}

sub _metric_header_as_key {
    return 'PAIR_ORIENTATION';
}

1;
