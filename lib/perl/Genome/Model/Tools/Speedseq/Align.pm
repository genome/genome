package Genome::Model::Tools::Speedseq::Align;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Speedseq::Align {
    is => 'Genome::Model::Tools::Speedseq::AlignBase',
    has_input => [
        fastq => {
            is => 'Text',
            doc => 'paired-end fastq file. if --paired flag is used then expected to be an interleaved paired-end fastq file, and fastq2 may be omitted. (can be gzipped)',
            example_values => ['in1.fq'],
            tool_bare_arg_position => '2',
            tool_input_file => 1,
        },
        fastq2 => {
            is => 'Text',
            doc => 'paired-end fastq file. (can be gzipped)',
            is_optional => 1,
            example_values => ['in2.fq'],
            tool_bare_arg_position => '3',
            tool_input_file => 1,
        }
    ],
    has_param => [
        read_group_header => {
            is => 'Text',
            doc => 'read group BAM header line',
            example_values => ['@RG\tID:id\tSM:samplename\tLB:lib'],
            tool_param_name => 'R',
        },
        paired => {
            is => 'Boolean',
            doc => 'first fastq file consists of interleaved paired-end sequences',
            is_optional => 1,
            tool_param_name => 'p',
        },
    ],
};

sub _tool_subcommand_name {
    return 'align';
}


#sub execute {
    #my $self = shift;
    #if ($self->paired) {
    #    if ($self->fastq2) {
    #        die('Expected single interleaved paired-end FASTQ!');
    #    }
    #    $options .= ' -p';
    #} else {
    #    push @input_files, $self->fastq2;
    #    $inputs_string .= ' '. $self->fastq2;
    #}

    # COMMAND
    # OUTPUTS
    #unless ($self->output_prefix) {
    #    $self->output_prefix($self->fastq);
    #}
#}

1;
