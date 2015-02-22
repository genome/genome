package Genome::InstrumentData::Command::Import::WorkFlow::SraToBam;

use strict;
use warnings;

use Genome;

require File::Basename;
require File::Spec;
use Try::Tiny;

class Genome::InstrumentData::Command::Import::WorkFlow::SraToBam { 
    is => 'Command::V2',
    has_input => [
        working_directory => {
            is => 'Text',
            doc => 'Destination directory for bam.',
        },
        sra_path => {
            is => 'Text',
            doc => 'Path of the SRA.',
        },
    ],
    has_output => {
        output_bam_path => {
            calculate_from => [qw/ working_directory sra_basename /],
            calculate => q| return File::Spec->join($working_directory, $sra_basename.'.bam'); |,
            doc => 'The path of the bam dumped from the SRA path.',
        },
        sra_basename => {
            calculate_from => [qw/ sra_path /],
            calculate => q| return File::Basename::basename($sra_path); |,
            doc => 'The basename of the SRA path.',
        },
    },
};

sub execute {
    my $self = shift;
    $self->debug_message('SRA to bam...');

    my $config_ok = $self->_check_ncbi_config;
    return unless $config_ok;

    my $dump_ok = $self->_dump_bam_from_sra;
    return if not $dump_ok;

    my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;
    my $flagstat = $helpers->validate_bam($self->output_bam_path);
    unless ($flagstat) {
        $self->error_message('Failed to validate bam.');
        return;
    }

    $self->debug_message('SRA to bam...done');
    return 1;
}

sub _check_ncbi_config {
    my $self = shift;
    $self->debug_message('Check for NCBI config file...');

    my $ncbi_config_file = $ENV{HOME}.'/.ncbi/user-settings.mkfg';
    if (-s $ncbi_config_file) {
        $self->debug_message('Check for NCBI config file...OK');
        return 1;
    }
    else {
        $self->error_message("No NCBI config file ($ncbi_config_file) found. "
            ."Please run 'perl /usr/bin/sra-configuration-assistant' to set it up. "
            ."This file is required for most NCBI SRA operations.");
        return;
    }
}

sub _dump_bam_from_sra {
    my $self = shift;
    $self->debug_message('Dump bam from SRA...');


    my $sra_path = $self->sra_path;
    $self->debug_message("SRA path: $sra_path");
    my $output_bam_path = $self->output_bam_path;
    $self->debug_message("Output bam path: $output_bam_path");

    # Capturing stderr of command in case the command fails so that an
    # informative error message can be recorded.  Using a temp file instead of
    # a temp var to protect against a potentially large output.
    my $stderr = $output_bam_path.'.err';

    try {
        Genome::Sys->shellcmd(
            cmd             => "/usr/bin/sam-dump --primary --unaligned $sra_path | samtools view -h -b -S -",
            redirect_stdout => $output_bam_path,
            redirect_stderr => $stderr,
        );
    }
    catch {
        $self->error_message('Caught exception from shellcmd: '. $_);
        my $err = Genome::Sys->open_file_for_reading($stderr);
        while (my $line = $err->getline) {
            $self->error_message("STDERR: $line");
        }
        return;
    };

    $self->debug_message('Dump bam from SRA...OK');
    return 1;
}

1;

