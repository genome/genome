package Genome::InstrumentData::Command::Import::WorkFlow::FastqsToBam;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
require File::Basename;
require List::Util;

class Genome::InstrumentData::Command::Import::WorkFlow::FastqsToBam { 
    is => 'Command::V2',
    has_input => [
        working_directory => {
            is => 'Text',
            doc => 'Detination directory for fastqs.',
        },
        fastq_paths => { 
            is => 'Text',
            is_many => 1,
            doc => 'Paths of the fastq[s] to convert to a bam.',
        },
        sample_name => {
            is => 'Text',
            doc => 'The sample name to use and then derive the read group name.',
        },
    ],
    has_output => [ 
        bam_path => {
            calculate_from => [qw/ working_directory sample_name /],
            calculate => q( return $working_directory.'/'.$sample_name.'.bam'; ),
            doc => 'The path of the bam.',
        },
    ],
};

sub execute {
    my $self = shift;
    $self->status_message('Fastqs to bam...');

    my $fastq_to_bam_ok = $self->_fastqs_to_bam;
    return if not $fastq_to_bam_ok;

    my $verfiy_read_count_ok = $self->_verify_read_counts;
    return if not $verfiy_read_count_ok;

    $self->status_message('Fastqs to bam...done');
    return 1;
}

sub _fastqs_to_bam {
    my $self = shift;
    $self->status_message('Run picard fastq to sam...');

    my $sample_name = $self->sample_name;
    my $read_group_name = $sample_name.'-extlibs';
    my @fastqs = $self->fastq_paths;
    $self->status_message("Fastq 1: $fastqs[0]");
    my $bam_path = $self->bam_path;
    my $cmd = "gmt picard fastq-to-sam --fastq $fastqs[0] --output $bam_path --quality-format Standard --sample-name $sample_name --read-group-name $read_group_name";
    if ( $fastqs[1] ) {
        $self->status_message("Fastq 2: $fastqs[1]") if $fastqs[1];
        $cmd .= ' --fastq2 '.$fastqs[1]
    }
    $self->status_message("Bam path: $bam_path");

    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv or not -s $bam_path ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to run picard fastq to sam!');
        return;
    }

    $self->status_message('Run picard fastq to sam...done');
    return 1;
}

sub _verify_read_counts {
    my $self = shift;
    $self->status_message('Verify read count...');

    my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;
    my $flagstat = $helpers->run_flagstat($self->bam_path);
    return if not $flagstat;
    $self->status_message('Bam read count: '.$flagstat->{total_reads});

    my $counts = $helpers->load_read_count_for_fastq_paths($self->fastq_paths);
    return if not $counts;
    my $fastq_read_count = List::Util::sum(values %$counts);
    $self->status_message('Fastq read count: '.$fastq_read_count);

    if ( $flagstat->{total_reads} != $fastq_read_count ) {
        $self->error_message('Fastq and bam rread counts do not match!');
        return;
    }

    $self->status_message('Verify read count...done');
    return 1;
}

1;

