package Genome::Model::Tools::Fastx::QualityStats;

use strict;
use warnings;

use Genome;
use Genome::Sys;
use File::Basename;

class Genome::Model::Tools::Fastx::QualityStats {
    is => ['Genome::Model::Tools::Fastx'],
    has_constant => {
        fastx_tool => { value => 'fastx_quality_stats' },
    },
    has_input => [
        fastq_file => {
            is => 'Text',
            doc => 'The fastq file to produce quality stats for',
        },
    ],
    has_output => [
        stats_file => {
            is => 'Text',
            doc => 'The output stats file',
            is_optional => 1,
        }
    ],
};

sub execute {
    my $self = shift;
    unless (Genome::Sys->validate_file_for_reading($self->fastq_file)) {
        $self->error_message('Failed to validate fastq file for read access '. $self->fastq_file .":  $!");
        die($self->error_message);
    }
    my @suffix = qw/fq fastq txt bam/;
    my ($basename,$dirname,$suffix) = File::Basename::fileparse($self->fastq_file,@suffix);
    $basename =~ s/\.$//;
    unless ($self->stats_file) {
        $self->stats_file($dirname .'/'. $basename .'.stats');
    }
    unless (Genome::Sys->validate_file_for_writing($self->stats_file)) {
        $self->error_message('Failed to validate stats file for write access '. $self->stats_file .":  $!");
        die($self->error_message);
    }
    #fastx_quality_stats (as of 0.0.7) won't process from tmp with -i and -o
    #my $cmd = $self->fastx_tool_path .' < '. $self->fastq_file .' > '. $self->stats_file;
    my $cmd = "gmt sx -input ".$self->fastq_file.":type=sanger -output file=-:type=illumina | ".$self->fastx_tool_path.' > '.$self->stats_file;
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$self->fastq_file],
        output_files => [$self->stats_file],
    );
    return 1;
}

1;
