package Genome::Model::Tools::SpeedSeq::AlignBase;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::SpeedSeq::AlignBase {
    is => 'Genome::Model::Tools::SpeedSeq::Base',
    is_abstract => 1,
    has_input => [
        reference_fasta => {
            is => 'Text',
            doc => 'fasta file (indexed with bwa)',
            example_values => ['reference.fa'],
            tool_bare_arg_position => '1',
            tool_input_file => 1,
        },
    ],
    has_param => [
        output_prefix => {
            is => 'Text',
            doc => 'output prefix',
            is_optional => 1,
            example_values => ['in1.fq','in.realign'],
            tool_param_name => 'o',
        },
        paired => {
            is => 'Boolean',
            doc => 'first fastq file consists of interleaved paired-end sequences',
            is_optional => 1,
            tool_param_name => 'p',
        },
        threads => {
            is => 'Integer',
            doc => 'threads',
            default_value => 1,
            is_optional => 1,
            tool_param_name => 't',
        },
        temp_directory => {
            is => 'Text',
            doc => 'temp directory',
            example_values => ['./output_prefix.XXXXXXXXXXXX'],
            is_optional => 1,
            tool_param_name => 'T',
        },
        insert_size => {
            is => 'Text',
            doc => 'specify the mean, standard deviation (10% of the mean if absent), max(4 sigma from the mean if absent) and min of the insert size distribution.  FR orientation only. Inferred by default.',
            is_optional => 1,
            example_values => ['FLOAT[,FLOAT[,INT[,INT]]]'],
            tool_param_name => 'I',
        },
        include_duplicates =>{
            is => 'Boolean',
            doc => 'include duplicates in splitters and discordants',
            is_optional => 1,
            tool_param_name => 'i',
        },
        maximum_split_alignments => {
            is => 'Integer',
            doc => 'maximum number of split alignments for a read to be included in splitter file',
            is_optional => 1,
            example_values => ['2'],
            tool_param_name => 'c',
        },
        minimum_non_overlapping_bp => {
            is => 'Integer',
            doc => 'minimum non-overlapping base pairs between two alignments for a read to be included in splitter file',
            is_optional => 1,
            example_values => ['20'],
            tool_param_name => 'm',
        },
        sort_memory => {
            is => 'Integer',
            doc => 'amount of memory in GB to be used for sorting',
            is_optional => 1,
            example_values => ['20'],
            tool_param_name => 'M',
        },
    ],
};

sub _output_extensions {
    my $class = shift;
    return qw/
                 bam
                 bam.bai
                 discordants.bam
                 discordants.bam.bai
                 splitters.bam
                 splitters.bam.bai
             /;
}

sub _resolve_output_files {
    my $self = shift;
    my @output_files = map { $self->output_prefix .'.'. $_ } $self->_output_extensions;
    $self->output_files(\@output_files);
    return @output_files;
}


1;
