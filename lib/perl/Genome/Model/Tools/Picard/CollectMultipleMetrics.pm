package Genome::Model::Tools::Picard::CollectMultipleMetrics;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Picard::CollectMultipleMetrics {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        input_file   => {
            is  => 'String',
            doc => 'The SAM/BAM files to run on.  File type is determined by suffix.',
        },
        output_basename  => {
            is  => 'String',
            doc => 'Basename to write output metrics files to.',
        },
        program_list => {
            is_optional => 1,
            doc => 'Comma-delimited ist of picard metrics collection programs to run.',
            default_value => join(',',
                'CollectAlignmentSummaryMetrics',
                'CollectInsertSizeMetrics',
                'QualityScoreDistribution',
                'MeanQualityByCycle'
                )
        },
        reference_sequence => {
            is_optional => 1,
        },
        stop_after => {
            is  => 'Integer',
            doc => 'Stop after processing N reads, mainly for debugging. Default value: 0.',
            is_optional   => 1,
        },
        assume_sorted => {
            is => 'Boolean',
            default_value => 1,
            is_optional => 1,
        },
    ],
};

sub help_brief {
    'Takes an input BAM and reference sequence and runs one or more Picard metrics modules at the same time to cut down on I/O. Currently all programs are run with default options and fixed output extesions, but this may become more flexible in future.';
}

sub help_detail {
    return <<EOS
    For Picard documentation of this command see:
    http://picard.sourceforge.net/command-line-overview.shtml#CollectMultipleMetrics
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

    my @args = $self->_translate_args(
        args => [qw(
            input_file
            output_basename
            stop_after
            reference_sequence
            assume_sorted
            )],

        translations => {
            input_file => 'INPUT',
            output_basename => 'OUTPUT',
        });

    my $programs = $self->program_list;
    $programs =~ s/ //g;
    my @programs = split(',', $self->program_list);
    push @args, map {sprintf "PROGRAM=%s", $_} @programs;

    return @args;
}

sub _validate_params {
    my $self = shift;
    $self->enforce_minimum_version('1.40');
}

sub execute {
    my $self = shift;
    return $self->basic_execute;
}

1;
