package Genome::Model::Tools::Picard::CollectMultipleMetrics;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Picard::CollectMultipleMetrics {
    is  => 'Genome::Model::Tools::Picard::Base',


    has_input => [
        input_file   => {
            is  => 'String',
            doc => 'The SAM/BAM files to run on.  File type is determined by suffix.',
            picard_param_name => 'INPUT',
        },
        output_basename  => {
            is  => 'String',
            doc => 'Basename to write output metrics files to.',
            picard_param_name => 'OUTPUT',
        },
        program_list => {
            is_optional => 1,
            is_many => 1,
            picard_param_name => 'PROGRAM',
            doc => 'List of picard metrics collection programs to run. Default: ' .
                join(',',
                    'CollectAlignmentSummaryMetrics',
                    'CollectInsertSizeMetrics',
                    'QualityScoreDistribution',
                    'MeanQualityByCycle'),
        },
        reference_sequence => {
            is_optional => 1,
            picard_param_name => 'REFERENCE_SEQUENCE',
        },
        stop_after => {
            is  => 'Integer',
            doc => 'Stop after processing N reads, mainly for debugging. Default value: 0.',
            is_optional   => 1,
            picard_param_name => 'STOP_AFTER',
        },
        assume_sorted => {
            is => 'Boolean',
            default_value => 1,
            is_optional => 1,
            picard_param_name => 'ASSUME_SORTED',
        },
    ],
};

sub help_brief {
    'Takes an input BAM and reference sequence and runs one or more Picard metrics modules at the same time to cut down on I/O. Currently all programs are run with default options and fixed output extesions, but this may become more flexible in future.';
}

sub help_detail {
    return <<EOS
    For Picard documentation of this command see:
    http://broadinstitute.github.io/picard/command-line-overview.html#CollectMultipleMetrics
EOS
}

sub _jar_name {
    return 'CollectMultipleMetrics.jar';
}

sub _java_class_name {
    return 'net.sf.picard.analysis.CollectMultipleMetrics';
}

sub _cmdline_args {
    my $self = shift;
    my @args = $self->SUPER::_cmdline_args;
    my @programs = $self->program_list;
    if (@programs) {
        # by default, adding PROGRAM=xyz /appends/ to the default list picard
        # already has. PROGRAM=null clears the default, so we prepend that to
        # the argument list.
        unshift @args, 'PROGRAM=null';
    }
    return @args;
}

sub _validate_params {
    my $self = shift;
    $self->enforce_minimum_version('1.40');
}

1;
