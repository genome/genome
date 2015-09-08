package Genome::Model::Tools::Picard::CollectGcBiasMetrics;

use strict;
use warnings;

use Genome;
use File::Basename;


class Genome::Model::Tools::Picard::CollectGcBiasMetrics {
    is  => 'Genome::Model::Tools::Picard::Base',
    has_input => [
        input_file   => {
            is  => 'String',
            doc => 'The BAM file to run on.',
            picard_param_name => 'INPUT',
        },
        output_file  => {
            is  => 'String',
            doc => 'The output metrics file',
            picard_param_name => 'OUTPUT',
        },
        refseq_file  => {
            is  => 'String',
            doc => 'The reference sequence file',
            picard_param_name => 'REFERENCE_SEQUENCE',
        },
        chart_output => {
            is  => 'String',
            doc => 'The PDF file to render the chart to. Default is GC_bias_chart.pdf in output_file dir',
            is_optional => 1,
            picard_param_name => 'CHART_OUTPUT',
            default_value => 'GC_bias_chart.pdf',
        },
        summary_output => {
            is  => 'String',
            doc => 'The text file to write summary metrics to. Default is GC_bias_summary.txt',
            is_optional => 1,
            picard_param_name => 'SUMMARY_OUTPUT',
            default_value => 'GC_bias_summary.txt',
        },
        window_size  => {
            is  => 'Integer',
            doc => 'The size of windows on the genome that are used to bin reads. Default value: 100',
            default_value => 100,
            is_optional   => 1,
            picard_param_name => 'WINDOW_SIZE',
        },
        min_genome_fraction => {
            is  => 'Number',
            doc => 'For summary metrics, exclude GC windows that include less than this fraction of the genome. Default value: 1.0E-5.',
            default_value => '1.0E-5',
            is_optional   => 1,
            picard_param_name => 'MINIMUM_GENOME_FRACTION',
        },
        max_records_in_ram => {
            is => 'Integer',
            doc => 'The number of alignment records to store in RAM before spilling to disk.',
            default_value => 500000,
            is_optional => 1,
            picard_param_name => 'MAX_RECORDS_IN_RAM',
        },
        assume_sorted => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'Assume that the BAM file is sorted, regardless of what the header says',
            picard_param_name => 'ASSUME_SORTED',
        },
        metric_accumulation_level => {
            is => 'Text',
            is_many => 1,
            is_optional => 1,
            picard_param_name => 'METRIC_ACCUMULATION_LEVEL',
        },
    ],
};

sub help_brief {
    'Tool to collect GC bias metrics from a BAM file.';
}

sub help_detail {
    return <<EOS
    For Picard documentation of this command see:
    http://broadinstitute.github.io/picard/command-line-overview.html#CollectGcBiasMetrics
EOS
}

sub _jar_name {
    return 'CollectGcBiasMetrics.jar';
}

sub _java_class {
    return qw(picard analysis CollectGcBiasMetrics);
}

1;
