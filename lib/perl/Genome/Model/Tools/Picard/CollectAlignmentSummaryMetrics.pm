package Genome::Model::Tools::Picard::CollectAlignmentSummaryMetrics;

use strict;
use warnings;

use Genome;
use Carp qw(confess);

class Genome::Model::Tools::Picard::CollectAlignmentSummaryMetrics {
    is  => 'Genome::Model::Tools::Picard::Base',
    has_input => [
        input_file => {
            is  => 'String',
            doc => 'The SAM/BAM files to run on.  File type is determined by suffix.',
            picard_param_name => 'INPUT',
        },
        output_file => {
            is  => 'String',
            doc => 'The output summary file',
            picard_param_name => 'OUTPUT',
        },
        refseq_file => {
            is  => 'String',
            doc => 'The reference sequence file',
            picard_param_name => 'REFERENCE_SEQUENCE',
        },
        assume_sorted => {
            is  => 'Boolean',
            doc => 'coordinate sorted beforehand or not, default 1',
            default_value => 1,
            is_optional   => 1,
            picard_param_name => 'ASSUME_SORTED',
        },
        max_insert_size => {
            is  => 'Integer',
            doc => 'Paired end reads above this insert size will be considered chimeric along with inter-chromosomal pairs. default 100000',
            default_value => 100000,
            is_optional   => 1,
            picard_param_name => 'MAX_INSERT_SIZE',
        },
        adapter_sequence => {
            is  => 'String',
            doc => 'adapter sequence to exclude',
            is_optional => 1,
            picard_param_name => 'ADAPTER_SEQUENCE',
        },
        is_bisulfite_sequenced => {
            is  => 'Boolean',
            doc => 'Whether the SAM or BAM file consists of bisulfite sequenced reads, default 0',
            default_value => 0,
            is_optional   => 1,
            picard_param_name => 'IS_BISULFITE_SEQUENCED',
        },
        metric_accumulation_level => {
            is => 'Text',
            is_many => 1,
            is_optional => 1,
            default_value => ['ALL_READS'],
            picard_param_name => 'METRIC_ACCUMULATION_LEVEL',
        }
    ],
};

sub help_brief {
    'Tool to collect alignment summary metrics from a SAM/BAM file.';
}

sub help_detail {
    return <<EOS
    For Picard documentation of this command see:
    http://broadinstitute.github.io/picard/command-line-overview.html#CollectAlignmentSummaryMetrics
EOS
}

sub _validate_params {
    my $self = shift;

    for my $type (qw(input refseq)) {
        my $property_name = $type.'_file';
        unless ($self->$property_name and -s $self->$property_name) {
            confess $self->error_message("$property_name is invalid");
        }
    }
}

sub _jar_name {
    return 'CollectAlignmentSummaryMetrics.jar';
}

sub _java_class {
    return qw(picard analysis CollectAlignmentSummaryMetrics);
}

sub _metric_header_as_key {
    return 'CATEGORY';
}

1;
