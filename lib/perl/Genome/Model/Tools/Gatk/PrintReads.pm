package Genome::Model::Tools::Gatk::PrintReads;

use strict;
use warnings;

use Genome;

require File::Basename;

class Genome::Model::Tools::Gatk::PrintReads {
    doc => "Run GATK with the 'PrintReads' tool",
    is => 'Genome::Model::Tools::Gatk::Base',
    has => [
        input_bams => {
            is_input => 1,
            is_many => 1,
            is => 'Text',
            doc => 'The path to the original bam(s) you would like to assess',
        },
        reference_fasta => {
            is_input => 1,
            is => 'Text',
            doc => "The path to the reference fasta you would like to run against" ,
        },
        output_bam => {
            is_input => 1,
            is_output => 1,
            is => 'Text',
            doc => 'The path to where you would like to create the output bam file',
        },
        bqsr => {
            is_input => 1,
            is => 'Text',
            is_optional => 1,
            doc => 'The input covariates table file which enables on-the-fly base quality score recalibration',
        },
        downsample_coverage => {
            is_input => 1,
            is_optional => 1,
            is => 'Number',
            doc => 'Downsample BAM to desired coverage',
        },
        number => {
            is_input => 1,
            is_optional => 1,
            is => 'Integer',
            doc => 'Print the first n reads from the file, discarding the rest',
        },
        platform => {
            is_input => 1,
            is_optional => 1,
            is => 'Text',
            doc => 'Exclude all reads with this platform from the output',
        },
        read_group => {
            is_input => 1,
            is_optional => 1,
            is => 'Text',
            doc => 'Exclude all reads with this read group from the output',
        },
        sample_file => {
            is_input => 1,
            is_many => 1,
            is_optional => 1,
            is => 'Text',
            doc => 'File containing a list of samples (one per line). Can be specified multiple times',
        },
        sample_name => {
            is_input => 1,
            is_many => 1,
            is_optional => 1,
            is => 'Text',
            doc => 'Sample name to be included in the analysis. Can be specified multiple times.',
        },
        simplify => {
            is_input => 1,
            is_optional => 1,
            is => 'Boolean',
            doc => 'Simplify all reads. Erase all extra attributes in the read but keep the read group information.',
        },
    ],
};

sub help_brief {
    "Run GATK with the 'PrintReads' tool"
}

sub help_synopsis {
    return <<EOS
    gmt gatk print-reads --input-bams input1.bam,input2.bam --reference-fasta ref.fasta --output-bam output.bam --bqsr recalibration_report.grp
EOS
}

sub execute {
    my $self = shift;

    unless ($self->_check_inputs) {
        return;
    }
    my $command = $self->print_reads_command;

    unless (Genome::Sys->shellcmd(cmd => $command)) {
        die $self->error_message("Failed to execute $command");
    }

    # Rename the bam index. It gets named without a .bam in the name
    my ($bam_basename, $bam_dirname) = File::Basename::fileparse($self->output_bam);
    $bam_basename =~ s/\.bam$//;
    Genome::Sys->rename($bam_dirname.'/'.$bam_basename.'.bai', $bam_dirname.'/'.$bam_basename.'.bam.bai');

    return 1;
}

sub print_reads_command {
    my $self = shift;
    my $gatk_command = $self->base_java_command;
    $gatk_command .= " -R " . $self->reference_fasta;
    $gatk_command .= " -T PrintReads";
    $gatk_command .= " -o " . $self->output_bam;
    $gatk_command .= " -I $_ " for $self->input_bams;

    $gatk_command .= " --BQSR " . $self->bqsr if $self->bqsr;

    $gatk_command .= " --downsample_coverage " . $self->downsample_coverage if $self->downsample_coverage;
    $gatk_command .= " --number " . $self->number if $self->number;
    $gatk_command .= " --platform " . $self->platform if $self->platform;
    $gatk_command .= " --readGroup " . $self->read_group if $self->read_group;
    $gatk_command .= " --sample_file " . $self->sample_file if $self->sample_file;
    $gatk_command .= " --sample_name " . $self->sample_name if $self->sample_name;
    $gatk_command .= " --simplify " . $self->simplify if $self->simplify;

    return $gatk_command;
}

sub _check_inputs {
    my $self = shift;

    Genome::Sys->validate_file_for_reading($_) for $self->input_bams;
    Genome::Sys->validate_file_for_reading($self->reference_fasta);
    Genome::Sys->validate_file_for_reading($self->bqsr) if $self->bqsr;
    Genome::Sys->validate_file_for_reading($_) for $self->sample_file;
    Genome::Sys->validate_file_for_writing($self->output_bam);

    return 1;
}

1;
