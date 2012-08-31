package Genome::Model::Tools::RefCov::WholeGenome;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::RefCov::WholeGenome {
    is => ['Genome::Model::Tools::RefCov'],
    has_input => [
        merged_stats_file => {
            is_optional => 0,
        },
        merge_by => {
            is_optional => 1,
            default_value => 'transcript',
        },
        evaluate_gc_content => {
            default_value => 1,
            is_optional => 1,
        },
        roi_normalized_coverage => {
            default_value => 1,
            is_optional => 1,
        },
        genome_normalized_coverage => {
            default_value => 1,
            is_optional => 1,
        },
        min_depth_filter => {
            default_value => '30',
            is_optional => 1,
        },
        print_headers => {
            default_value => 1,
            is_optional => 1,
        },

    ],
};

sub help_brief {
    "The default for running RefCov on WholeGenome BAMs with 30x depth.  Includes merged results by transcript and normalized coverage by genome and ROI.",
}

sub execute {
    my $self = shift;
    unless ($self->print_roi_coverage) {
        die('Failed to print ROI coverage!');
    }
    return 1;
}

1;  # end of package
