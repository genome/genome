package Genome::Model::Tools::RefCov::RnaSeq;

use strict;
use warnings;

use Genome;

my $DEFAULT_MAXIMUM_ALIGNMENT_COUNT = 100_000;

class Genome::Model::Tools::RefCov::RnaSeq {
    is => ['Genome::Model::Tools::RefCov'],
    has => [
        print_headers => {
            is_optional => 1,
            default_value => 1,
        },
        merged_stats_file => {
            is_optional => 0,
        },
        merge_by => {
            is_optional => 1,
            default_value => 'transcript',
        },
        # This was becoming too expensive to run on large BAM files - jrw 09/04/2012
        #alignment_count => {
        #    default_value => 1,
        #    is_optional => 1,
        #},
        print_min_max => {
            default_value => 1,
            is_optional => 1,
        },
        maximum_alignment_count => {
            default_value => $DEFAULT_MAXIMUM_ALIGNMENT_COUNT,
            is_optional => 1,
        }
    ],
};

sub help_brief {
    "Default settings customized for running ref-cov on typical RNA-seq BAMs.  Additional RNA-seq outputs over the standard mode include an alignment count with minimum and maximum depth as well as transcript-level metrics given properly formatted BED files(see gtf-to-bed).",
}

sub execute {
    my $self = shift;
    unless ($self->print_roi_coverage) {
        die('Failed to print ROI coverage!');
    }
    return 1;
}

1;
