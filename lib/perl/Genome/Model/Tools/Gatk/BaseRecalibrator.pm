package Genome::Model::Tools::Gatk::BaseRecalibrator;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Gatk::BaseRecalibrator {
    doc => "Run GATK with the 'BaseRecalibrator' tool",
    is => [qw/ Genome::Model::Tools::Gatk::Base Genome::Model::Tools::Gatk::WithNumberOfCpuThreads /],
    has => [
        input_bam => {
            is_input => 1,
            is => 'Text',
            doc => 'The path to the original bam you would like to assess',
        },
        reference_fasta => {
            is_input => 1,
            is => 'Text',
            doc => "The path to the reference fasta you would like to run against" ,
        },
        known_sites => {
            is_input => 1,
            is_many => 1,
            is_optional => 1,
            is => 'Text',
            doc => 'A database of known polymorphic sites to skip over in the recalibration algorithm',
        },
        output_recalibration_table => {
            is_input => 1,
            is_output => 1,
            is => 'Text',
            doc => 'The path to where you would like to create the output recalibration table file',
        },
    ],
};

sub help_brief {
    "Run GATK with the 'BaseRecalibrator' tool"
}

sub help_synopsis {
    return <<EOS
    gmt gatk base-recalibator --input-bam my_reads.bam --reference-fasta my.fa --known-sites path/to/sites.vcf,path/to/more/sites.vcf --output-recalibration-table recal_data.grp
EOS
}

sub execute {
    my $self = shift;

    unless ($self->_check_inputs) {
        return;
    }
    my $command = $self->base_recalibrator_command;

    unless (Genome::Sys->shellcmd(cmd => $command)) {
        die $self->error_message("Failed to execute $command");
    }

    return 1;
}

sub base_recalibrator_command {
    my $self = shift;

    my $gatk_command = $self->base_java_command;
    $gatk_command .= " -T BaseRecalibrator";
    $gatk_command .= " -I " . $self->input_bam;
    $gatk_command .= " -R " . $self->reference_fasta;
    $gatk_command .= " -knownSites $_" for $self->known_sites;
    $gatk_command .= " -o " . $self->output_recalibration_table;
    $gatk_command .= $self->number_of_cpu_threads_param_for_java_command;

    return $gatk_command;
}

sub _check_inputs {
    my $self = shift;

    Genome::Sys->validate_file_for_reading($self->input_bam);
    Genome::Sys->validate_file_for_reading($self->reference_fasta);
    Genome::Sys->validate_file_for_reading($_) for $self->known_sites;
    Genome::Sys->validate_file_for_writing($self->output_recalibration_table);

    return 1;
}

1;
