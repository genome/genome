package Genome::Model::Tools::Picard::CollectAlignmentSummaryMetrics;

use strict;
use warnings;

use Genome;


class Genome::Model::Tools::Picard::CollectAlignmentSummaryMetrics {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        input_file => {
            is  => 'String',
            doc => 'The SAM/BAM files to run on.  File type is determined by suffix.',
        },
        output_file => {
            is  => 'String',
            doc => 'The output summary file',
        },
        refseq_file => {
            is  => 'String',
            doc => 'The reference sequence file',
        },
        assume_sorted => {
            is  => 'Boolean',
            doc => 'coordinate sorted beforehand or not, default 1',
            default_value => 1,
            is_optional   => 1,
        },
        max_insert_size => {
            is  => 'Integer',
            doc => 'Paired end reads above this insert size will be considered chimeric along with inter-chromosomal pairs. default 100000',
            default_value => 100000,
            is_optional   => 1,
        },
        adapter_sequence => {
            is  => 'String',
            doc => 'adapter sequence to exclude',
            is_optional => 1,
        },
        is_bisulfite_sequenced => {
            is  => 'Boolean',
            doc => 'Whether the SAM or BAM file consists of bisulfite sequenced reads, default 0',
            default_value => 0,
            is_optional   => 1,
        },
    ],
};

sub help_brief {
    'Tool to collect alignment summary metrics from a SAM/BAM file.';
}

sub help_detail {
    return <<EOS
    For Picard documentation of this command see:
    http://picard.sourceforge.net/command-line-overview.shtml#CollectAlignmentSummaryMetrics
EOS
}

sub execute {
    my $self = shift;

    for my $type (qw(input refseq)) {
        my $property_name = $type.'_file';
        unless ($self->$property_name and -s $self->$property_name) {
            $self->error_message("$property_name is invalid");
            return;
        }
    }

    my $cmd = $self->picard_path .'/CollectAlignmentSummaryMetrics.jar net.sf.picard.analysis.CollectAlignmentSummaryMetrics';
    $cmd   .= ' OUTPUT='. $self->output_file  .' INPUT='. $self->input_file .' REFERENCE_SEQUENCE='. $self->refseq_file;
    
    if ($self->assume_sorted) {
        $cmd .= ' ASSUME_SORTED=true';
    }
    if ($self->max_insert_size) {
        $cmd .= ' MAX_INSERT_SIZE='. $self->max_insert_size;
    }
    if ($self->adapter_sequence) {
        $cmd .= ' ADAPTER_SEQUENCE='. $self->adapter_sequence;
    }
    if ($self->is_bisulfite_sequenced) {
        $cmd .= ' IS_BISULFITE_SEQUENCED=true';
    }
    $self->run_java_vm(
        cmd          => $cmd,
        input_files  => [$self->input_file],
        output_files => [$self->output_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}

sub _metric_header_as_key {
    return 'CATEGORY';
}

1;
