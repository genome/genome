package Genome::Model::Tools::SpeedSeq::Align;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::SpeedSeq::Align {
    is => 'Genome::Model::Tools::SpeedSeq::AlignBase',
    has_input => [
        fastq => {
            is => 'Text',
            doc => 'paired-end fastq file. if -p flag is used then expected to be an interleaved paired-end fastq file, and fastq2 may be omitted. (can be gzipped)',
            example_values => ['in1.fq'],
        },
        fastq2 => {
            is => 'Text',
            doc => 'paired-end fastq file. (can be gzipped)',
            is_optional => 1,
            example_values => ['in2.fq'],
        }
    ],
    has_param => [
        read_group_header => {
            is => 'Text',
            doc => 'read group BAM header line',
            example_values => ['@RG\tID:id\tSM:samplename\tLB:lib'],
        },
    ],
};


sub execute {
    my $self = shift;

    # OPTIONS aka params
    my $options = $self->_resolve_options_string('-R "'. $self->read_group_header .'"');

    # INPUTS
    my @input_files;
    my $inputs_string = $self->reference_fasta .' '. $self->fastq;
    push @input_files, $self->reference_fasta;
    push @input_files, $self->fastq;
    if ($self->config_file) {
        push @input_files, $self->config_file;
    }
    if ($self->paired) {
        if ($self->fastq2) {
            die('Expected single interleaved paired-end FASTQ!');
        }
        $options .= ' -p';
    } else {
        push @input_files, $self->fastq2;
        $inputs_string .= ' '. $self->fastq2;
    }

    # COMMAND
    my $cmd = $self->speedseq_path .' align '. $options .' '. $inputs_string;

    # OUTPUTS
    my $basename = $self->fastq;
    if ($self->output_prefix) {
        $basename = $self->output_prefix;
    }
    my @output_files = $self->_resolve_output_file_paths($basename);
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => \@input_files,
        output_files => \@output_files,
        skip_if_output_is_present => 0,
        keep_dbh_connection_open => 0,
    );
    $self->output_files(\@output_files);

    return 1;
}

1;
