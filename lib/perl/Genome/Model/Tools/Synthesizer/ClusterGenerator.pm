package Genome::Model::Tools::Synthesizer::ClusterGenerator;

use strict;
use warnings;

use Genome;

my $DEFAULT_ZENITH = '5';
my $DEFAULT_MIN_DEPTH = '1';

class Genome::Model::Tools::Synthesizer::ClusterGenerator {
    is => ['Genome::Model::Tools::Synthesizer::Base'],
    has_input => [
        sized_bam_file => {
        	is => 'Text',
            doc => 'Input file of size-selected and filtered alignments in BAM format',
        },
        zenith_depth => {
            is => 'Text',
            is_output=> 1,
            doc => 'A minimum zenith(maximum depth) size to retain a cluster',
            default_value => $DEFAULT_ZENITH,
        },
        minimum_depth => {
            is => 'Text',
            is_output=> 1,
            doc => 'Minimum depth to filter coverage',
            default_value => $DEFAULT_MIN_DEPTH,
        },
        stats_file => {
            is => 'Text',
            is_output=> 1,
            doc => 'Calculate statistics across clusters and print to file ',
        },
        bed_file => {
            is => 'Text',
            is_output=> 1,
            doc => 'The output BED format file of clusters',
        },
    ],
};

sub help_brief {
"Run the Synthesizer Cluster-Generator module to identify regions of contiguous coverage called \"Clusters\" across the genome and identify putative sncRNA species ";
}

sub help_detail {
"Run the Synthesizer Cluster-Generator module to identify regions of contiguous coverage called \"Clusters\" across the genome and identify putative sncRNA species";
}


sub execute {
    my $self = shift;

    my $cmd = Genome::Model::Tools::BioSamtools::ClusterCoverage->create(
        bam_file => $self->sized_bam_file,
        minimum_zenith => $self->zenith_depth,
        minimum_depth => $self->minimum_depth,
        stats_file => $self->stats_file,
        bed_file => $self->bed_file,
    );
    unless ($cmd->execute) {
        die $self->error_message('Failed to run cluster-coverage command.');
    }

    return 1;
}

1;

