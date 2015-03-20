package Genome::InstrumentData::Command::Import::WorkFlow::FastqsToBam;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
require File::Basename;
require List::Util;
use Try::Tiny;

class Genome::InstrumentData::Command::Import::WorkFlow::FastqsToBam { 
    is => 'Command::V2',
    has_input => [
        working_directory => {
            is => 'Text',
            doc => 'Destination directory for fastqs.',
        },
        fastq_paths => { 
            is => 'Text',
            is_many => 1,
            doc => 'Paths of the fastq[s] to convert to a bam.',
        },
        sample_name => {
            is => 'Text',
            doc => 'The sample name to use in the RG header.',
        },
        library_name => {
            is => 'Text',
            doc => 'The library name to use as the IDs in the SAM RG and LB headers.',
        },
    ],
    has_output => [ 
        output_bam_path => {
            is => 'Text',
            calculate_from => [qw/ working_directory sample_name /],
            calculate => q( return $working_directory.'/'.$sample_name.'.bam'; ),
            doc => 'The path of the bam.',
        },
    ],
};

sub execute {
    my $self = shift;
    $self->debug_message('Fastqs to bam...');

    my $unarchive_if_necessary = $self->_unarchive_fastqs_if_necessary;
    return if not $unarchive_if_necessary;

    my $fastq_to_bam_ok = $self->_fastqs_to_bam;
    return if not $fastq_to_bam_ok;

    my $verify_bam_ok = $self->_verify_bam;
    return if not $verify_bam_ok;

    my $cleanup_ok = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->remove_paths_and_auxiliary_files($self->fastq_paths);
    return if not $cleanup_ok;

    $self->debug_message('Fastqs to bam...done');
    return 1;
}

sub _unarchive_fastqs_if_necessary {
    my $self = shift;
    $self->debug_message('Unarchive fastqs if necessary...');

    my @new_fastq_paths;
    for my $fastq_path ( $self->fastq_paths ) {
        if ( $fastq_path !~ /\.gz$/ ) {
            $self->debug_message('Unarchive not necessary for '.$fastq_path);
            push @new_fastq_paths, $fastq_path;
            next;
        }
        $self->debug_message('Unarchiving: '.$fastq_path);
        
        my $success = try {
            Genome::Sys->shellcmd(cmd => "gunzip $fastq_path");
        }
        catch {
            $self->error_message($_) if $_;
            $self->error_message('Failed to gunzip fastq!');
            return;
        };
        return if not $success;

        my $unarchived_fastq_path = $fastq_path;
        $unarchived_fastq_path =~ s/\.gz$//;
        $self->debug_message("Unarchived fastq: $unarchived_fastq_path");
        if ( not -s $unarchived_fastq_path ) {
            $self->error_message('Unarchived fastq does not exist!');
            return;
        }
        push @new_fastq_paths, $unarchived_fastq_path;
        unlink $fastq_path;
    }
    $self->fastq_paths(\@new_fastq_paths);

    $self->debug_message('Unarchive fastqs if necessary...');
    return 1;
}

sub _fastqs_to_bam {
    my $self = shift;
    $self->debug_message('Run picard fastq to sam...');

    my @fastqs = $self->fastq_paths;
    $self->debug_message("Fastq 1: $fastqs[0]");
    my $output_bam_path = $self->output_bam_path;
    my %fastq_to_sam_params = (
        fastq => $fastqs[0],
        output => $output_bam_path,
        quality_format => 'Standard',
        sample_name => $self->sample_name,
        library_name => $self->library_name,
        read_group_name => $self->library_name,
        use_version => '1.113',
    );
    if ( $fastqs[1] ) {
        $self->debug_message("Fastq 2: $fastqs[1]");
        $fastq_to_sam_params{fastq2} = $fastqs[1];
    }
    $self->debug_message("Bam path: $output_bam_path");

    my $cmd = Genome::Model::Tools::Picard::FastqToSam->create(%fastq_to_sam_params);
    if ( not $cmd ) {
        $self->error_message('Failed to create sam to fastq command!');
        return;
    }
    my $success = try {
        $cmd->execute;
    }
    catch {
        $self->error_message($_) if $_;
        $self->error_message('Failed to run picard fastq to sam!');
        return;
    };
    return if not $success;

    if ( not -s $output_bam_path ) {
        $self->error_message('Ran picard fastq to sam, but bam path does not exist!');
        return;
    }

    $self->debug_message('Run picard fastq to sam...done');
    return 1;
}

sub _verify_bam {
    my $self = shift;
    $self->debug_message('Verify bam...');

    my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;

    my $flagstat = $helpers->validate_bam($self->output_bam_path);
    return if not $flagstat;

    $self->debug_message('Bam read count: '.$flagstat->{total_reads});

    $self->debug_message('Verify bam...done');
    return 1;
}

1;

