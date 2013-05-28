package Genome::Model::Tools::Gatk::RealignerTargetCreator;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Gatk::RealignerTargetCreator {
    doc => "Run GATK with the 'RealignerTargetCreator' tool",
    is => 'Genome::Model::Tools::Gatk',
    has_input => [
        known => {
            is => 'Text',
            doc => 'The file of known indels',
            is_optional => 1,
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
        max_interval_size => {
            is => 'Integer',
            doc => "maximum interval size; any intervals larger than this value will be dropped",
            example_values => [500],
        },
        min_reads_at_locus => {
            is => 'Integer',
            doc => "minimum reads at a locus to enable using the entropy calculation",
            example_values => [4],
        },
        mismatch_fraction => {
            is => 'Number',
            doc => 'fraction of base qualities needing to mismatch for a position to have high entropy',
            example_values => [0.0],
        },
        window_size => {
            is => 'Integer',
            doc => 'window size for calculating entropy or SNP clusters',
            example_values => [10],
        },
    ],
};

sub help_brief {
    "Run GATK with the 'RealignerTargetCreator' tool"
}

sub help_synopsis {
    return <<EOS
    gmt gatk realignr-indels-target-creator --known indels.vcf --input-bam my_existing.bam --reference-fasta my.fa --max-interval-size 500 --min-reads-at-locus 4 --mismatch-fraction 0.0 --window-size 10 --output-intervals out.intervals
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
    if ($self->known) {
       $gatk_command .= " --known " . $self->known;
    }
    $gatk_command .= " --maxIntervalSize " . $self->max_interval_size;
    $gatk_command .= " --minReadsAtLocus " . $self->min_reads_at_locus;
    $gatk_command .= " --mismatchFraction " . $self->mismatch_fraction;
    $gatk_command .= " --windowSize " . $self->window_size;
    $gatk_command .= " -I " . $self->input_bam;
    $gatk_command .= " -R " . $self->reference_fasta;
    $gatk_command .= " -o ". $self->output_intervals;
    return $gatk_command;
}

sub _check_inputs {
    my $self = shift;

    if ($self->known) {
        Genome::Sys->validate_file_for_reading($self->known);
    }
    Genome::Sys->validate_file_for_reading($self->input_bam);
    Genome::Sys->validate_file_for_reading($self->reference_fasta);
    Genome::Sys->validate_file_for_writing($self->output_intervals);

    return 1;
}

1;
