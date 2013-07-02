package Genome::InstrumentData::Command::Import::WorkFlow::ConvertFastqsToBam;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
require File::Basename;
require List::Util;
require List::MoreUtils;

class Genome::InstrumentData::Command::Import::WorkFlow::ConvertFastqsToBam { 
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
    $self->status_message('Convert fastqs to bam...');

    my $convert_ok = $self->_convert_fastqs_to_bam;
    return if not $convert_ok;

    #my $flagstat = Genome::InstrumentData::Command::Helpers->run_flagstat();

    $self->status_message('Convert fastqs to bam...done');
    return 1;
}

sub _convert_fastqs_to_bam {
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

1;

