package Genome::InstrumentData::Command::Import::Workflow::SraToBam;

use strict;
use warnings;

use Genome;

require Cwd;
require File::Basename;

class Genome::InstrumentData::Command::Import::Workflow::SraToBam { 
    is => 'Command::V2',
    has_input => [
        working_directory => {
            is => 'Text',
            doc => 'Destination directory for bam.',
        },
        sra_path => {
            is => 'Text',
            is_many => 1,
            doc => 'Path of the SRA.',
        },
    ],
    has_output => [
        bam_path => {
            calculate_from => [qw/ working_directory /],
            calculate => q( return $working_directory.'/dumped_from_sra.bam'; ),
            doc => 'The path of the bam dumped from the SRA path.',
        },
    ],
};

sub execute {
    my $self = shift;
    $self->status_message('Dump bam from SRA...');

    my $dump_ok = $self->_dump_bam_from_sra;
    return if not $dump_ok;

    my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;
    my $flagstat = $helpers->run_flagstat($self->sorted_bam_path);
    return if not $flagstat;

    $self->status_message('Dump bam from SRA...done');
    return 1;
}

sub _dump_bam_from_sra {
    my $self = shift;

    $self->status_message('Check for NCBI config file...');
    my $ncbi_config_file = $ENV{HOME}.'/.ncbi/user-settings.mkfg';
    if ( not -s $ncbi_config_file ) {
        $self->error_message("No NCBI config file ($ncbi_config_file) found. Please run 'perl /usr/bin/sra-configuration-assistant' to set it up. This file is required for most NCBI SRA operations.");
        return
    }
    $self->status_message('Check for NCBI config file...done');

    my $sra_path = $self->sra_path;
    $self->status_message('SRA path: '.$sra_path);

    $self->status_message('Check SRA database...');
    my $dbcc_file = $sra_path.'.dbcc';
    $self->status_message('DBCC file: '.$dbcc_file);
    my $cwd = Cwd::getcwd();
    my ($source_sra_basename, $source_sra_directory) = File::Basename::fileparse($sra_path);
    chdir $source_sra_directory;
    my $cmd = "/usr/bin/sra-dbcc $source_sra_basename &> $dbcc_file";
    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv or not -s $dbcc_file ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to run sra dbcc!');
        return;
    }
    my @dbcc_lines = eval{ Genome::Sys->read_file($dbcc_file); };
    if ( not @dbcc_lines ) {
        $self->error_message('Failed to read SRA DBCC file! ');
        return;
    }
    my $sra_has_primary_alignment_info = grep { $_ =~ /PRIMARY_ALIGNMENT/ } @dbcc_lines;
    chdir $cwd;
    $self->status_message('Check SRA database...done');
    
    $self->status_message('Dump aligned bam...');
    my $aligned_bam = $self->_tmp_dir.'/aligned.bam';
    my $unsorted_bam = $aligned_bam; # set now, override if there is a unaligned bam
    $self->status_message('Aligned bam: '.$aligned_bam);
    $cmd = "/usr/bin/sam-dump --primary $sra_path | samtools view -h -b -S - > $aligned_bam";
    $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv or not -s $aligned_bam ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to run sra sam dump aligned bam!');
        return;
    }
    $self->status_message('Dump aligned bam...done');

    if ( $sra_has_primary_alignment_info ) { # unaligned are already dumped above if there is no alignment info
        $self->status_message('Dump unaligned from sra to fastq...');
        my $unaligned_fastq = $self->_tmp_dir.'/unaligned.fastq';
        $self->status_message("Unaligned fastq: $unaligned_fastq");
        $cmd = "/usr/bin/fastq-dump --unaligned --origfmt $sra_path --stdout > $unaligned_fastq";
        $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
        if ( not $rv ) {
            $self->error_message($@) if $@;
            $self->error_message('Failed to run sra sam dump unaligned fastq!');
            return;
        }
        $self->status_message('Dump unaligned from sra to fastq...done');

        if ( -s $unaligned_fastq ) {
            $self->status_message('Convert unaligned fastq to bam...');
            my $unaligned_bam = $unaligned_fastq.'.bam';
            my $cmd = "gmt picard fastq-to-sam --fastq $unaligned_fastq --output $unaligned_bam --quality-format Standard --sample-name ".$self->sample->name;
            my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
            if ( not $rv or not -s $unaligned_bam ) {
                $self->error_message($@) if $@;
                $self->error_message('Failed to run sam fastq to sam on unaligned fastq!');
                return;
            }
            $self->status_message('Convert unaligned fastq to bam...done');
            unlink($unaligned_fastq);

            $self->status_message('Add bam from unaligned fastq to unsorted bam...');
            $unsorted_bam = $self->_tmp_dir.'/unsorted.bam';
            $cmd = "samtools merge $unsorted_bam $aligned_bam $unaligned_bam";
            $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
            if ( not $rv ) {
                $self->error_message($@) if $@;
                $self->error_message('Failed to run samtools view!');
                return;
            }
            $self->status_message('Add bam from unaligned fastq to unsorted bam...done');
            unlink($unaligned_bam);
        }
    }

    $self->status_message('Transfer SRA file...done');
    return 1;
}

1;

