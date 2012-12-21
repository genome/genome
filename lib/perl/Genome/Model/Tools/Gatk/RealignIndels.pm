package Genome::Model::Tools::Gatk::RealignIndels;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Gatk::RealignIndels {
    doc => "Run GATK with the 'IndelRealigner' tool",
    is => 'Genome::Model::Tools::Gatk',
    has => [
        target_intervals => {
            is => 'Text',
            doc => 'The file of indels around which you wish to do realignment',
        },
        output_realigned_bam => {
            is => 'Text',
            doc => 'The path to where you would like the realigned output bam',
        },
        input_bam => {
            is => 'Text',
            doc => 'The path to the original bam you would like to be realigned',
        },
        target_intervals_are_sorted => {
            is => 'Boolean',
            doc => 'If set to false, pass along --targetIntervalsAreNotSorted',
            default => 0,
        },
        reference_fasta => {
            is => 'Text',
            doc => "Reference Fasta" ,
        },
        index_bam => {
            is => 'Boolean',
            default => 1,
            doc => 'Index the bam after alignment.'
        }
    ],
};

sub help_brief {
    "Run GATK with the 'IndelRealigner' tool"
}

sub help_synopsis {
    return <<EOS
    gmt gatk realign-indels --target-intervals some.bed --output-realigned-bam my_output_realigned.bam --input-bam my_existing.bam --reference-fasta my.fa
EOS
}

sub execute {
    my $self = shift;

    $self->_check_inputs;
    my $command = $self->indel_realigner_command;

    unless (Genome::Sys->shellcmd(cmd => $command)) {
        die $self->error_message("Failed to execute $command");
    }

    if ($self->index_bam) {
        my $aligned_bam = $self->output_realigned_bam;
        die $self->error_message("Couldn't find realigned bam at $aligned_bam!") unless -f $aligned_bam;

        my $rv = Genome::Model::Tools::Sam::IndexBam->execute(bam_file => $aligned_bam);
        die $self->error_message("Failed to run gmt sam index-bam on $aligned_bam") unless $rv->result == 1;
    }

    return 1;
}

sub indel_realigner_command {
    my $self = shift;
    my $gatk_command = $self->base_java_command;
    $gatk_command .= " -T IndelRealigner";
    $gatk_command .= " -targetIntervals " . $self->target_intervals;
    $gatk_command .= " -o " . $self->output_realigned_bam;
    $gatk_command .= " -I " . $self->input_bam;
    $gatk_command .= " -R " . $self->reference_fasta;
    unless ($self->target_intervals_are_sorted) {
        $gatk_command .= " --targetIntervalsAreNotSorted";
    }
    return $gatk_command;
}

sub _check_inputs {
    my $self = shift;

    Genome::Sys->validate_file_for_reading($self->target_intervals);
    Genome::Sys->validate_file_for_reading($self->input_bam);
    Genome::Sys->validate_file_for_reading($self->reference_fasta);
    Genome::Sys->validate_file_for_writing($self->output_realigned_bam);

    return 1;
}

1;
