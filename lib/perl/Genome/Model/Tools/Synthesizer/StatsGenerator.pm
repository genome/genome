package Genome::Model::Tools::Synthesizer::StatsGenerator;

use strict;
use warnings;

use Genome;

my $DEFAULT_CUTOFF = '2';

class Genome::Model::Tools::Synthesizer::StatsGenerator {
    is => 'Genome::Model::Tools::Synthesizer::Base',
    has_input => [
        coverage_stats_file => {
            is_output => 1,
            is        => 'Text',
            doc       => 'Input stats file from ClusterCoverage',
        },
        sized_bam_file => {
            is        => 'Text',
            doc       => 'Input size-specific BAM file of alignments.Make sure the BAM index .bai file is also in the same directory. To index, run \'samtools index\'',
            is_output => 1
        },
        head_bin_flagstat_file => {
            is        => 'Text',
            doc       => 'Input flagstat file of the normalization bin',
            is_output => 1
        },
        output_stats_file => {
            is        => 'Text',
            is_output => 1,
            doc       => 'Output TSV file containing alignment statistics for the clusters ',
        },
        output_clusters_file => {
            is        => 'Text',
            is_output => 1,
            doc       => 'Output BED file containing coordinates of clusters in BED format (sorted by depth) ',
        },
        output_subclusters_file => {
            is        => 'Text',
            is_output => 1,
            doc       => 'Output BED file of "Subclusters" for each Cluster in the input BED file',
        },
        output_subcluster_intersect_file => {
            is        => 'Text',
            is_output => 1,
            doc       => 'Output TSV file of Subclusters that map with existing clusters',
        },
        subcluster_min_mapzero => {
            is            => 'Text',
            is_output     => 1,
            doc           => 'Minimum %MapZero Alignments to call subclusters',
            default_value => $DEFAULT_CUTOFF,
        },
    ],


};

sub help_brief {
    "Run the Synthesizer Stats-Generator module to calculate alignment as well coverage statistics for each Cluster. ";
}

sub help_detail {
    "Run the Synthesizer Stats-Generator module to calculate alignment as well coverage statistics for each Cluster.The output is a TSV file containing statistics info. This module also outputs regions where there equally best alignments for a cluster ; also called  \"sub-clusters\" ";
}

sub output_file_headers {
    return (
        "Cluster", "Chr", "Start", "Stop", "Avg Depth", "Zenith Depth",
        "Length of Raw Cluster", "# Positive Strand", "# Negative Strand",
        "Log Normalization -head bin", "Log Normalization -per bin",
        "% Mismatches", "ZeroMM", "1MM", "2MM", "3MM", "4MM",
        "% 1st Pos MM ", "Avg MapQ", "Std Dev Map Q", "%Zero MapQ",
        "Avg BaseQ", "Major Subcluster Loci",
    );
}

sub resolve_bam_file { shift->sized_bam_file }
sub resolve_input_flagstat_file { shift->head_bin_flagstat_file }

sub execute {
    my $self = shift;

    Genome::Model::SmallRna::Command::StatsGenerator->class
        || die "Can't load Genome::Model::SmallRna::Command::StatsGenerator";
    return $self->Genome::Model::SmallRna::Command::StatsGenerator::_execute_body();
}

1;
