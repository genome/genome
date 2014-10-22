package Genome::InstrumentData::Command::Import::WorkFlow::SraToBam;

use strict;
use warnings;

use Genome;
use Try::Tiny;
use File::Copy qw(move);

require Cwd;
require File::Basename;

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
        library => {
            is => 'Genome::Library',
            doc => 'The library name to use and then derive the read group name.',
        },
    ],
    has_output => [
        output_bam_path => {
            calculate_from => [qw/ sra_path /],
            calculate => q( return $sra_path.'.bam'; ),
            doc => 'The path of the bam dumped from the SRA path.',
        },
    ],
};

sub execute {
    my $self = shift;

    return try {
        $self->debug_message('Check for NCBI config file...');
        my $config_ok = $self->check_ncbi_config;
        return unless $config_ok;

        $self->debug_message('Dump bam from SRA...');
        my $dump_ok = $self->_dump_bam_from_sra;
        return if not $dump_ok;

        my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;
        my $flagstat = $helpers->validate_bam($self->output_bam_path);
        unless ($flagstat) {
            $self->error_message('Failed to validate bam.');
            return;
        }

        $self->debug_message('Dump bam from SRA...done');
        return 1;
    }
    catch {
        $self->error_message('Dumping BAM from SRA failed: '.$_);
        return;
    };
}


sub _dump_bam_from_sra {
    my $self = shift;
    
    my $aligned_bam = $self->working_directory.'/aligned.bam';
    my $dump_aligned_bam_ok = $self->dump_aligned_bam($self->sra_path, $aligned_bam);
    if ( $dump_aligned_bam_ok and (-s $aligned_bam) ) {
        $self->debug_message('Dump aligned bam...done');
    }
    else {
        $self->error_message('Failed to run sra sam dump aligned bam!');
        return;
    }

    $self->debug_message('Check SRA database...');
    my $sra_has_primary_alignment_info = $self->check_sra_database;
    $self->debug_message('Check SRA database...done');

    $self->debug_message('Dump aligned bam...');
    if ( $sra_has_primary_alignment_info ) {
        # if primary alignment info exists, only aligned are dumped above.
        return $self->add_unaligned_reads_to_bam($aligned_bam);
    }
    else {
        unless ( move($aligned_bam, $self->output_bam_path) ) {
            $self->error_message(
                'Failed to move aligned bam to output bam path. '
                .'Source: %s, Destination %s, Error Status %s',
                $aligned_bam, $self->output_bam_path, $!);
        }
    }

    return 1;
}

sub check_ncbi_config {
    my $self = shift;

    $self->debug_message('Check for NCBI config file...done');
    my $ncbi_config_file = $ENV{HOME}.'/.ncbi/user-settings.mkfg';

    if (-s $ncbi_config_file) {
        return 1;
    }
    else {
        $self->error_message("No NCBI config file ($ncbi_config_file) found. "
            ."Please run 'perl /usr/bin/sra-configuration-assistant' to set it up. "
            ."This file is required for most NCBI SRA operations.");
        return;
    }
}

sub run_sra_dbcc {
    my $self = shift;
    my ($source_sra_basename, $source_sra_directory, $dbcc_file) = @_;

    my $cwd = Cwd::getcwd();
    chdir($source_sra_directory) or die "Failed to chdir('$source_sra_directory')";
    my $cmd = "/usr/bin/sra-dbcc $source_sra_basename 2>&1";
    my $sra_dbcc_ok = $self->do_shellcmd_with_stdout($cmd, $dbcc_file);
    chdir($cwd) or die "Failed to chdir('$cwd')";
    return $sra_dbcc_ok;
}

sub check_sra_database {
    my $self = shift;
    my $sra_path = $self->sra_path;
    my $dbcc_file = $sra_path.'.dbcc';

    $self->debug_message('DBCC file: '.$dbcc_file);
    my ($source_sra_basename, $source_sra_directory) = File::Basename::fileparse($sra_path);
    my $run_ok = $self->run_sra_dbcc($source_sra_basename, $source_sra_directory, $dbcc_file);
    unless ($run_ok) {
        $self->error_message('Failed to run sra-dbcc');
        return;
    }

    my @dbcc_lines = Genome::Sys->read_file($dbcc_file);
    my $sra_has_primary_alignment_info = grep { $_ =~ /PRIMARY_ALIGNMENT/ } @dbcc_lines;
    return $sra_has_primary_alignment_info;
}


sub add_unaligned_reads_to_bam {
    my $self = shift;
    my ($aligned_bam) = @_;

    $self->debug_message('Dump unaligned from sra to fastq...');
    my $unaligned_fastq = $self->working_directory.'/unaligned.fastq';
    if ( $self->dump_unaligned_fastq($self->sra_path, $unaligned_fastq) ) {
        $self->debug_message('Dump unaligned from sra to fastq...done');
    }
    else {
        $self->error_message('Failed to run sra fastq-dump !');
        return;
    }

    if ( -s $unaligned_fastq ) {
        return unless $self->merge_unaligned_fastq_into_bam($unaligned_fastq, $aligned_bam);
    }
    else {
        unless ( move($aligned_bam, $self->output_bam_path) ) {
            $self->error_message(
                'Failed to move aligned bam to output bam path. '
                .'Source: %s, Destination %s, Error Status %s',
                $aligned_bam, $self->output_bam_path, $!);
        }
    }

    unlink($unaligned_fastq);
}

sub merge_unaligned_fastq_into_bam {
    my $self = shift;
    my ($unaligned_fastq, $aligned_bam) = @_;

    $self->debug_message('Convert unaligned fastq to bam...');
    my $unaligned_bam = $unaligned_fastq.'.bam';

    my $conversion_ok = $self->convert_fastq_to_bam(
        $self->library->sample->name,
        $unaligned_fastq, $unaligned_bam);
    if ($conversion_ok) {
        $self->debug_message('Convert unaligned fastq to bam...done');
    }
    else {
        $self->error_message('Failed to convert unaligned fastq to bam.');
        return;
    }

    $self->debug_message('Add bam from unaligned fastq to unsorted bam...');
    my $merge_ok = $self->merge_bams($aligned_bam, $unaligned_bam, $self->output_bam_path);
    if ($merge_ok) {
        $self->debug_message('Add bam from unaligned fastq to unsorted bam...done');
    }
    else {
        $self->error_message('Failed to add bam from unaligned fastq to unsorted bam');
        return;
    }

    unlink($unaligned_bam);

    return 1;
}

sub do_shellcmd {
    my $self = shift;
    my ($cmd) = @_;

    my ($stdout, $stderr) = Genome::Sys->create_temp_file_path for qw(1 2);

    # Capturing output of command in case the command fails so that an
    # informative error message can be recorded.  Using temp files instead of
    # temp vars to protect against potentially large outputs.

    return try {
        Genome::Sys->shellcmd(
            cmd             => $cmd,
            redirect_stdout => $stdout,
            redirect_stderr => $stderr);
        return 1;
    }
    catch {
        $self->error_message('Caught exception from shellcmd: '. $_);

        my $out = Genome::Sys->open_file_for_reading($stdout);
        $self->error_message('STDOUT: '. $_) while $out->getline;

        my $err = Genome::Sys->open_file_for_reading($stderr);
        $self->error_message('STDERR: '. $_) while $err->getline;

        return;
    };
}

sub do_shellcmd_with_stdout {
    my $self = shift;
    my ($cmd, $output_path) = @_;

    my $stderr = join('.', $output_path, 'err');

    # Capturing stderr of command in case the command fails so that an
    # informative error message can be recorded.  Using a temp file instead of
    # a temp var to protect against a potentially large output.

    return try {
        Genome::Sys->shellcmd(
            cmd             => $cmd,
            redirect_stdout => $output_path,
            redirect_stderr => $stderr);
        return 1;
    }
    catch {
        $self->error_message('Caught exception from shellcmd: '. $_);
        my $err = Genome::Sys->open_file_for_reading($stderr);
        $self->error_message('STDERR: '. $_) while $err->getline;
        $self->error_message($_) while $err->getline;
        return;
    };
}

sub dump_aligned_bam {
    my $self = shift;
    my ($sra_path, $aligned_bam) = @_;

    my $command = "/usr/bin/sam-dump --primary $sra_path | samtools view -h -b -S -";
    return $self->do_shellcmd_with_stdout( $command, $aligned_bam);
}

sub dump_unaligned_fastq {
    my $self = shift;
    my ($sra_path, $unaligned_fastq) = @_;

    my $command = "/usr/bin/fastq-dump --unaligned --origfmt --stdout $sra_path";
    return $self->do_shellcmd_with_stdout( $command, $unaligned_fastq);
}

sub convert_fastq_to_bam {
    my $self = shift;
    my ($sample_name, $unaligned_fastq, $unaligned_bam) = @_;

    return try {
        my $fastq_to_sam = Genome::Model::Tools::Picard::FastqToSam->create(
            fastq           => "$unaligned_fastq",
            output          => "$unaligned_bam",
            sample_name     => "$sample_name",
            quality_format  => 'Standard');
        return $fastq_to_sam->execute;
    }
    catch {
        $self->error_message('Caught exception from '
            .'Genome::Model::Tools::Picard::FastqToSam: '. $_);
        return;
    };
}

sub merge_bams {
    my $self = shift;
    my ($a,$b,$dest) = @_;

    my $command = "samtools merge $dest $a $b";
    return $self->do_shellcmd($command);
}

1;
