package Genome::Model::Tools::Picard::CollectInsertSizeMetrics;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Picard::CollectInsertSizeMetrics {
    is => 'Genome::Model::Tools::Picard::Base',
    has_input => [
        input_file => {
            is => 'String',
            doc => 'The SAM/BAM files to run on.  File type is determined by suffix.',
            picard_param_name => 'INPUT',
        },
        output_file => {
            is => 'String',
            doc => 'File to write insert size metrics to.',
            picard_param_name => 'OUTPUT',
        },
        histogram_file => {
            # FIXME: we should make this optional
            is => 'String',
            doc => 'File to write insert size histogram chart to.',
            picard_param_name => 'HISTOGRAM_FILE',
        },
        tail_limit => {
            is => 'Integer',
            doc => 'When calculating mean and stdev stop when the bins in the tail of the distribution contain fewer than mode/TAIL_LIMIT items. This also limits how much data goes into each data category of the histogram. Default value: 10000.',
            is_optional => 1,
            picard_param_name => 'TAIL_LIMIT',
        },
        histogram_width => {
            is => 'Integer',
            doc => 'Explicitly sets the histogram width, overriding the TAIL_LIMIT option. Also, when calculating mean and stdev, only bins <= HISTOGRAM_WIDTH will be included.Default value: null.',
            is_optional => 1,
            picard_param_name => 'HISTOGRAM_WIDTH',
        },
        minimum_pct  => {
            is => 'Integer',
            doc => 'When generating the histogram, discard any data categories (out of FR, TANDEM, RF) that have fewer than this percentage of overall reads. (Range: 0 to 1) Default value: 0.01.',
            is_optional   => 1,
            picard_param_name => 'MINIMUM_PCT',
        },
        stop_after => {
            is => 'Integer',
            doc => 'Stop after processing N reads, mainly for debugging. Default value: 0.',
            is_optional   => 1,
            picard_param_name => 'STOP_AFTER',
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

sub _jar_name {
    return 'CollectInsertSizeMetrics.jar';
}

sub _java_class_name {
    return 'net.sf.picard.analysis.CollectInsertSizeMetrics';
}

sub _metric_header_as_key {
    return 'PAIR_ORIENTATION';
}

1;
