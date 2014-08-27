package Genome::Model::Tools::Gatk::RealignerTargetCreator;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Gatk::RealignerTargetCreator {
    doc => "Run GATK with the 'RealignerTargetCreator' tool",
    is => [qw/ Genome::Model::Tools::Gatk::Base Genome::Model::Tools::Gatk::WithNumberOfThreads /],
    has_input => [
        known => {
            is => 'Text',
            doc => 'The file of known indels',
            is_optional => 1,
            is_many => 1,
        },
        input_bam => {
            is => 'Text',
            doc => 'The path to the original bam you would like to be realigned',
        },
        reference_fasta => {
            is => 'Text',
            doc => "Reference Fasta" ,
        },
        output_intervals => {
            is => 'Text',
            doc => "File of intervals to target for realignment",
        },
    ],
};

sub help_brief {
    "Run GATK with the 'RealignerTargetCreator' tool"
}

sub help_synopsis {
    return <<EOS
    gmt gatk realigner-indels-target-creator --known indels.vcf --input-bam my_existing.bam --reference-fasta my.fa --output-intervals out.intervals
EOS
}

sub execute {
    my $self = shift;

    unless ($self->_check_inputs) {
        return;
    }
    my $command = $self->realigner_creator_command;

    unless (Genome::Sys->shellcmd(cmd => $command)) {
        die $self->error_message("Failed to execute $command");
    }

    return 1;
}

sub realigner_creator_command {
    my $self = shift;
    my $gatk_command = $self->base_java_command;
    $gatk_command .= " -T RealignerTargetCreator";
    my @known = $self->known;
    if (@known) {
       $gatk_command .= " --known ".join(" --known ", @known);
    }
    $gatk_command .= " -I " . $self->input_bam;
    $gatk_command .= " -R " . $self->reference_fasta;
    $gatk_command .= " -o ". $self->output_intervals;
    $gatk_command .= $self->number_of_threads_param_for_java_command;
    return $gatk_command;
}

sub _check_inputs {
    my $self = shift;

    my @known = $self->known;
    if (@known) {
        for my $k (@known) {
            Genome::Sys->validate_file_for_reading($k);
        }
    }
    Genome::Sys->validate_file_for_reading($self->input_bam);
    Genome::Sys->validate_file_for_reading($self->reference_fasta);
    Genome::Sys->validate_file_for_writing($self->output_intervals);

    return 1;
}

1;
