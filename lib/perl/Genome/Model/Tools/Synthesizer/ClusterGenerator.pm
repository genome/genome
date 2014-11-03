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


    my $cmd = 'genome-perl5.10 -S gmt bio-samtools cluster-coverage --bam-file='. $self->sized_bam_file .' --minimum-zenith='. $self->zenith_depth .' --minimum-depth='. $self->minimum_depth .' --stats-file='. $self->stats_file .' --bed-file='. $self->bed_file ;
   # needs to be changed to G::M::T::RefCov::ClusterCoverage and released with RefCov
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$self->sized_bam_file],
        output_files => [$self->bed_file,$self->stats_file],
        skip_if_output_is_present => 0,
    );

    return 1;
}

1;

__END__

