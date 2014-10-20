package Genome::InstrumentData::Command::Import::WorkFlow::SraToBam;

use strict;
use warnings;

use Genome;
use Try::Tiny;
use IO::File;
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
        $self->debug_message('Dump bam from SRA...');
        my $dump_ok = $self->_dump_bam_from_sra;
        return if not $dump_ok;

        my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;
        my $flagstat = $helpers->validate_bam($self->output_bam_path);
        return if not $flagstat;

        $self->debug_message('Dump bam from SRA...done');
        return 1;
    }
    catch {
        $self->debug_message('Dumping BAM from SRA failed: '.$_);
        return;
    };
}

sub _dump_bam_from_sra {
    my $self = shift;

    $self->debug_message('Check for NCBI config file...');
    my $ncbi_config_file = $ENV{HOME}.'/.ncbi/user-settings.mkfg';
    if ( not -s $ncbi_config_file ) {
        $self->error_message("No NCBI config file ($ncbi_config_file) found. Please run 'perl /usr/bin/sra-configuration-assistant' to set it up. This file is required for most NCBI SRA operations.");
        return
    }
    $self->debug_message('Check for NCBI config file...done');

    my $sra_path = $self->sra_path;
    $self->debug_message('SRA path: '.$sra_path);

    $self->debug_message('Check SRA database...');
    my $dbcc_file = $sra_path.'.dbcc';
    $self->debug_message('DBCC file: '.$dbcc_file);
    my $cwd = Cwd::getcwd();
    my ($source_sra_basename, $source_sra_directory) = File::Basename::fileparse($sra_path);
    chdir($source_sra_directory) or die "Failed to chdir('$source_sra_directory')";
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
    chdir($cwd) or die "Failed to chdir('$cwd')";
    $self->debug_message('Check SRA database...done');
    
    $self->debug_message('Dump aligned bam...');
    my $aligned_bam = ( $sra_has_primary_alignment_info )
        ? $self->working_directory.'/aligned.bam'
        : $self->output_bam_path;

    my $dump_aligned_bam_ok = $self->dump_aligned_bam($sra_path, $aligned_bam);
    if ( $dump_aligned_bam_ok and (-s $aligned_bam) ) {
        $self->debug_message('Dump aligned bam...done');
    }
    else {
        $self->error_message('Failed to run sra sam dump aligned bam!');
        return;
    }

    if ( $sra_has_primary_alignment_info ) {
        # if primary alignment info exists, only aligned are dumped above.
        $self->debug_message('Dump unaligned from sra to fastq...');

        my $unaligned_fastq = $self->working_directory.'/unaligned.fastq';
        if ( $self->dump_unaligned_fastq($sra_path, $unaligned_fastq) ) {
            $self->debug_message('Dump unaligned from sra to fastq...done');
        }
        else {
            $self->error_message('Failed to run sra fastq-dump !');
            return;
        }

        if ( -s $unaligned_fastq ) {
            my $unaligned_bam = $unaligned_fastq.'.bam';

            $self->debug_message('Convert unaligned fastq to bam...');
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
            my $bam_path = $self->output_bam_path;
            $cmd = "samtools merge $bam_path $aligned_bam $unaligned_bam";
            $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
            if ( not $rv ) {
                $self->error_message($@) if $@;
                $self->error_message('Failed to run samtools view!');
                return;
            }
            $self->debug_message('Add bam from unaligned fastq to unsorted bam...done');
            unlink($unaligned_bam);
        }

        else {
            unless ( move($aligned_bam, $self->output_bam_path) ) {
                $self->error_message( sprintf(
                    'Failed to move aligned bam to output bam path. '
                    .'Source: %s, Destination %s, Error Status %s',
                    $aligned_bam, $self->output_bam_path, $!));
            }
        }

        unlink($unaligned_fastq);
    }

    return 1;
}

sub do_shellcmd {
    my $self = shift;
    my ($cmd) = @_;

    my ($stdout, $stderr) = Genome::Sys->create_temp_file_path for qw(1 2);

    return try {
        Genome::Sys->shellcmd(
            cmd             => $cmd,
            redirect_stdout => $stdout,
            redirect_stderr => $stderr);
        return 1;
    }
    catch {
        $self->error_message('Caught exception from shellcmd: '. $_);

        my $out = IO::File->new;
        $out->open($stdout, '<');
        $self->debug_message('STDOUT'. $_) while $out->getline;

        my $err = IO::File->new;
        $err->open($stderr, '<');
        $self->debug_message('STDERR'. $_) while $err->getline;

        return;
    };
}

sub do_shellcmd_with_stdout {
    my $self = shift;
    my ($cmd, $output_path) = @_;

    my $stderr = join('.', $output_path, 'err');

    return try {
        Genome::Sys->shellcmd(
            cmd             => $cmd,
            redirect_stdout => $output_path,
            redirect_stderr => $stderr);
        return 1;
    }
    catch {
        $self->error_message('Caught exception from shellcmd: '. $_);
        my $fh = IO::File->new;
        $fh->open($stderr, '<');
        $self->debug_message($_) while $fh->getline;
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

    my $command = "gmt picard fastq-to-sam "
        ."--fastq $unaligned_fastq "
        ."--output $unaligned_bam "
        ."--quality-format Standard --sample-name $sample_name";
    return $self->do_shellcmd($command);
}


1;
