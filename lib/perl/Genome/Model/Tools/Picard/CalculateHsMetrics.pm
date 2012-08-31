package Genome::Model::Tools::Picard::CalculateHsMetrics;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Picard::CalculateHsMetrics {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        input_file   => {
            is  => 'String',
            doc => 'The SAM/BAM files to run on.  File type is determined by suffix.',
        },
        output_file  => {
            is  => 'String',
            doc => 'The output metrics file',
        },
        bait_intervals  => {
            is  => 'String',
            doc => 'An interval list file that contains the locations of the baits used.',
        },
        target_intervals => {
            is  => 'String',
            doc => 'An interval list file that contains the locations of the targets',
        },
        reference_sequence => {
            is => 'String',
            doc => 'The reference sequence aligned to. Default value: null.',
            is_optional => 1,
        },
        per_target_coverage_file => {
            is => 'String',
            doc => '',
            is_optional => 1,
        },
        bait_set_name => {
            is => 'String',
            doc => 'Bait set name. If not provided it is inferred from the filename of the bait intervals.',
            is_optional => 1,
        },
        metric_accumulation_level => {
            is => 'String',
            doc => 'The level(s) at which to accumulate metrics. Possible values: {ALL_READS, SAMPLE, LIBRARY, READ_GROUP} This option may be specified 0 or more times. This option can be set to \'null\' to clear the default list.',
            is_optional => 1,
        },
    ],
};

sub help_brief {
    'Calculates a set of Hybrid Selection specific metrics from an aligned SAM or BAM file.';
}

sub help_detail {
    return <<EOS
    For Picard documentation of this command see:
    http://picard.sourceforge.net/command-line-overview.shtml#CalculateHsMetrics
EOS
}

sub execute {
    my $self = shift;

    my $cmd = $self->picard_path .'/CalculateHsMetrics.jar net.sf.picard.analysis.directed.CalculateHsMetrics OUTPUT='. $self->output_file
        .' INPUT='. $self->input_file .' BAIT_INTERVALS='. $self->bait_intervals .' TARGET_INTERVALS='. $self->target_intervals;
    
    my @input_files = ($self->input_file,$self->bait_intervals,$self->target_intervals);
    my @output_files = ($self->output_file);
    
    if ($self->reference_sequence) {
        push @input_files, $self->reference_sequence;
        $cmd .= ' REFERENCE_SEQUENCE='. $self->reference_sequence;
    }
    if ($self->per_target_coverage_file) {
        push @output_files, $self->per_target_coverage_file;
        $cmd .= ' PER_TARGET_COVERAGE='. $self->per_target_coverage_file;
    }
    if ($self->metric_accumulation_level) {
        $cmd .= ' METRIC_ACCUMULATION_LEVEL='. $self->metric_accumulation_level;
    }
    if ($self->bait_set_name) {
        if ($self->use_version > '1.52') {
            $cmd .= ' BAIT_SET_NAME='. $self->bait_set_name;
        } else {
            $self->warning_message('bait_set_name is not compatible with Picard version '. $self->use_version);
        }
    }

    $self->run_java_vm(
        cmd          => $cmd,
        input_files  => \@input_files,
        output_files => \@output_files,
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
