package Genome::Model::Tools::Speedseq::AlignBase;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Speedseq::AlignBase {
    is => 'Genome::Model::Tools::Speedseq::Base',
    is_abstract => 1,
    has_input => [
        reference_fasta => {
            is => 'Text',
            doc => 'fasta file (indexed with bwa)',
            example_values => ['reference.fa'],
        },
    ],
    has_param => [
        output_prefix => {
            is => 'Text',
            doc => 'output prefix',
            is_optional => 1,
            example_values => ['in1.fq'],
        },
        paired => {
            is => 'Boolean',
            doc => 'first fastq file consists of interleaved paired-end sequences',
            is_optional => 1,
        },
        threads => {
            is => 'Integer',
            doc => 'threads',
            default_value => 1,
            is_optional => 1,
        },
        temp_directory => {
            is => 'Text',
            doc => 'temp directory',
            example_values => ['./output_prefix.XXXXXXXXXXXX'],
            is_optional => 1,
        },
        insert_size => {
            is => 'Text',
            doc => 'specify the mean, standard deviation (10% of the mean if absent), max(4 sigma from the mean if absent) and min of the insert size distribution.  FR orientation only. Inferred by default.',
            is_optional => 1,
            example_values => ['FLOAT[,FLOAT[,INT[,INT]]]'],
        },
        include_duplicates =>{
            is => 'Boolean',
            doc => 'include duplicates in splitters and discordants',
            is_optional => 1,
        },
        maximum_split_alignments => {
            is => 'Integer',
            doc => 'maximum number of split alignments for a read to be included in splitter file',
            is_optional => 1,
            example_values => ['2'],
        },
        minimum_non_overlapping_bp => {
            is => 'Integer',
            doc => 'minimum non-overlapping base pairs between two alignments for a read to be included in splitter file',
            is_optional => 1,
            example_values => ['20'],
        },
        sort_memory => {
            is => 'Integer',
            doc => 'amount of memory in GB to be used for sorting',
            is_optional => 1,
            example_values => ['20'],
        },
    ],
    has_optional_output => [
        output_files => {
            is_many => 1,
        },
    ],
};

sub _output_extensions {
    my $class = shift;
    return qw/
                 .bam
                 .bam.bai
                 .discordants.bam
                 .discordants.bam.bai
                 .splitters.bam
                 .splitters.bam.bai
             /;
}

sub _resolve_options_string {
    my $self = shift;
    my $params = shift;

    unless ($params) {
        $params = '';
    }
    if ($self->output_prefix) {
        $params .= ' -o '. $self->output_prefix;
    }
    if ($self->threads) {
        $params .= ' -t '. $self->threads;
    }
    if ($self->temp_directory) {
        $params .= ' -T '. $self->temp_directory;
    }
    if ($self->insert_size) {
        $params .= ' -I '. $self->insert_size;
    }
    if ($self->include_duplicates) {
        $params .= ' -i ';
    }
    if ($self->maximum_split_alignments) {
        $params .= ' -c '. $self->maximum_split_alignments;
    }
    if ($self->minimum_non_overlapping_bp) {
        $params .= ' -m '. $self->minimum_non_overlapping_bp;
    }
    if ($self->sort_memory) {
        $params .= ' -M '. $self->sort_memory;
    }
    if ($self->config_file) {
        $params .= ' -K '. $self->config_file;
    }
    if ($self->verbose) {
        $params .= ' -v ';
    }
    return $params;
}

sub _resolve_output_file_paths {
    my $self = shift;
    my $basename = shift;
    my @output_files;
    for my $extension ($self->_output_extensions) {
        push @output_files, $basename . $extension;
    }
    return @output_files;
}

1;
