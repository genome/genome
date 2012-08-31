package Genome::Model::Tools::BioSamtools::ReadLengthDistribution;

use strict;
use warnings;

use Statistics::Descriptive;
use Genome;

class Genome::Model::Tools::BioSamtools::ReadLengthDistribution {
    is => 'Genome::Model::Tools::BioSamtools',
    has_input => [
        bam_file => {
            is => 'Text',
            doc => 'A path to a BAM format file of aligned capture reads',
        },
        summary_file => {
            is => 'Text',
            doc => 'A summary of the distribution metrics.',
        },
        read_length_histogram => {
            is => 'Text',
            doc => 'The output read length distribution file.',
        },
        alignment_length_histogram => {
            is => 'Text',
            doc => 'The output alignment length distribution file.',
        }
    ],
};

sub execute {
    my $self = shift;

    my @length_headers = ('Length','Count','Proportion');
    my $read_length_writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $self->read_length_histogram,
        separator => "\t",
        headers => \@length_headers,
    );
    unless ($read_length_writer) { die; }
    
    my $alignment_length_writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $self->alignment_length_histogram,
        separator => "\t",
        headers => \@length_headers,
    );
    unless ($alignment_length_writer) { die; }

    my @summary_headers = ('Type','Count','Minimum','Maximum','Mean','Standard_Deviation');
    my $summary_writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $self->summary_file,
        separator => "\t",
        headers => \@summary_headers,
    );
    unless ($summary_writer) { die; }
    
    my $refcov_bam  = Genome::Model::Tools::RefCov::Bam->create(bam_file => $self->bam_file );
    unless ($refcov_bam) {
        die('Failed to load bam file '. $self->bam_file);
    }
    
    my $bam  = $refcov_bam->bio_db_bam;
    my $header = $bam->header();

    my %read_stats;
    my $read_stats = Statistics::Descriptive::Sparse->new();
    my %alignment_stats;
    my $alignment_stats = Statistics::Descriptive::Sparse->new();

    while (my $align = $bam->read1()) {
        my $flag = $align->flag;
        my $read_length = $align->l_qseq;
        $read_stats->add_data($read_length);
        $read_stats{$read_length}++;
        unless ($flag & 4) {
            my $align_length = $align->calend - $align->pos;
            $alignment_stats{$align_length}++;
            $alignment_stats->add_data($align_length);
        }
    }

    # The total count is actually alignments, NOT reads...
    my %read_summary = (
        'Type' => 'Reads',
        'Count' => $read_stats->count,
        'Minimum' => $read_stats->min,
        'Maximum' => $read_stats->max,
        'Mean' => $read_stats->mean,
        'Standard_Deviation' => $read_stats->standard_deviation,
    );
    $summary_writer->write_one(\%read_summary);
    for my $key (sort {$a <=> $b} keys %read_stats) {
        my %data = (
            'Length' => $key,
            'Count' => $read_stats{$key},
            'Proportion' => ($read_stats{$key} / $read_stats->count),
        );
        $read_length_writer->write_one(\%data);
    }

    my %alignment_summary = (
        'Type' => 'Alignments',
        'Count' => $alignment_stats->count,
        'Minimum' => $alignment_stats->min,
        'Maximum' => $alignment_stats->max,
        'Mean' => $alignment_stats->mean,
        'Standard_Deviation' => $alignment_stats->standard_deviation,
    );
    $summary_writer->write_one(\%alignment_summary);
    for my $key (sort {$a <=> $b} keys %alignment_stats) {
        my %data = (
            'Length' => $key,
            'Count' => $alignment_stats{$key},
            'Proportion' => ($alignment_stats{$key} / $alignment_stats->count),
        );
        $alignment_length_writer->write_one(\%data);
    }
    return 1;
}


1;
