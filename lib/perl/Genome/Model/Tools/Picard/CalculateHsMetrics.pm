package Genome::Model::Tools::Picard::CalculateHsMetrics;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Picard::CalculateHsMetrics {
    is  => 'Genome::Model::Tools::Picard::Base',
    has_input => [
        input_file   => {
            is  => 'String',
            doc => 'The SAM/BAM files to run on.  File type is determined by suffix.',
            picard_param_name => 'INPUT',
        },
        output_file  => {
            is  => 'String',
            doc => 'The output metrics file',
            picard_param_name => 'OUTPUT',
        },
        bait_intervals  => {
            is  => 'String',
            doc => 'An interval list file that contains the locations of the baits used.',
            picard_param_name => 'BAIT_INTERVALS',
        },
        target_intervals => {
            is  => 'String',
            doc => 'An interval list file that contains the locations of the targets',
            picard_param_name => 'TARGET_INTERVALS',
        },
        reference_sequence => {
            is => 'String',
            doc => 'The reference sequence aligned to. Default value: null.',
            is_optional => 1,
            picard_param_name => 'REFERENCE_SEQUENCE',
        },
        per_target_coverage_file => {
            is => 'String',
            doc => '',
            is_optional => 1,
            picard_param_name => 'PER_TARGET_COVERAGE',
        },
        bait_set_name => {
            is => 'String',
            doc => 'Bait set name. If not provided it is inferred from the filename of the bait intervals.',
            is_optional => 1,
            picard_param_name => 'BAIT_SET_NAME',
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
    'Calculates a set of Hybrid Selection specific metrics from an aligned SAM or BAM file.';
}

sub help_detail {
    return <<EOS
    For Picard documentation of this command see:
    http://broadinstitute.github.io/picard/command-line-overview.html#CalculateHsMetrics
EOS
}

sub _jar_name {
    return 'CalculateHsMetrics.jar';
}

sub _java_class {
    return qw(picard analysis directed CalculateHsMetrics);
}

sub _shellcmd_extra_params {
    my $self = shift;

    my @input_files = grep {defined $_} ($self->reference_sequence, $self->input_file);
    my @output_files = grep {defined $_} ($self->per_target_coverage_file);

    return (
        input_files => \@input_files,
        output_files => \@output_files,
        skip_if_output_is_present => 0,
        );
}

sub _validate_params {
    my $self = shift;

    if ($self->bait_set_name && $self->version_at_most('1.52')) {
        $self->warning_message(
            sprintf 'bait_set_name is not compatible with Picard version %s, ignoring.',
                $self->use_version
            );

        $self->bait_set_name(undef);
    }
}

1;
