package Genome::Model::Tools::Picard::EstimateLibraryComplexity;

use strict;
use warnings;

use Genome;


class Genome::Model::Tools::Picard::EstimateLibraryComplexity {
    is  => 'Genome::Model::Tools::Picard::Base',
    has_input => [
        input_file => {
            is  => 'String',
            doc => 'One or more SAM/BAM files to combine and estimate library complexity from. Reads can be mapped or unmapped.',
            is_many => 1,
            picard_param_name => 'INPUT',
        },
        output_file => {
            is  => 'String',
            doc => 'The output metrics file from picard..',
            picard_param_name => 'OUTPUT',
        },
        min_identical_bases => {
            doc => 'The minimum number of bases at the starts of reads that must be identical for reads to be grouped together for duplicate detection. In effect total_reads / 4^max_id_bases reads will be compared at a time, so lower numbers will produce more accurate results but consume exponentially more memory and CPU.',
            default_value => 5,
            is_optional => 1,
            picard_param_name =>' MIN_IDENTICAL_BASES',
        },
        max_diff_rate => {
            doc => 'The maximum rate of differences between two reads to call them identical.',
            is_optional => 1,
            default_value => '0.03',
            picard_param_name => 'MAX_DIFF_RATE',
        },
        min_mean_quality => {
            doc => 'The minimum mean quality of the bases in a read pair for the read to be analyzed. Reads with lower average quality are filtered out and not considered in any calculations.',
            default_value => 20,
            is_optional => 1,
            picard_param_name => 'MIN_MEAN_QUALITY',
        },
        read_name_regex => {
            is => 'Text',
            doc => 'Regular expression that can be used to parse read names in the incoming SAM file. Read names are parsed to extract three variables: tile/region, x coordinate and y coordinate. These values are used to estimate the rate of optical duplication in order to give a more accurate estimated library size. The regular expression should contain three capture groups for the three variables, in order.',
            is_optional => 1,
            #default_value => '[a-zA-Z0-9\-\_]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*',
            picard_param_name => 'READ_NAME_REGEX',
        },
        optical_duplicate_pixel_distance => {
            doc => 'The maximum offset between two duplicte clusters in order to consider them optical duplicates. This should usually be set to some fairly small number (e.g. 5-10 pixels) unless using later versions of the Illumina pipeline that multiply pixel values by 10, in which case 50-100 is more normal.',
            default_value => 100,
            is_optional => 1,
            picard_param_name => 'OPTICAL_DUPLICATE_PIXEL_DISTANCE',
        },
    ],
};

sub help_brief {
    'Tool to estimate library complexity of raw reads(even unaligned) from a SAM/BAM file.';
}

sub help_detail {
    return <<EOS
    For Picard documentation of this command see:
    http://broadinstitute.github.io/picard/command-line-overview.html#EstimateLibraryComplexity
EOS
}

sub _jar_name {
    return 'EstimateLibraryComplexity.jar';
}

sub _java_class {
    return qw(picard sam EstimateLibraryComplexity);
}

sub _shellcmd_extra_params {
    my $self = shift;
    return (
        input_files => [$self->input_file],
        output_files => [$self->output_file],
        skip_if_output_is_present => 0,
        );
}

1;
