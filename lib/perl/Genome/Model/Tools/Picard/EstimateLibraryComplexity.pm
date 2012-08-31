package Genome::Model::Tools::Picard::EstimateLibraryComplexity;

use strict;
use warnings;

use Genome;


class Genome::Model::Tools::Picard::EstimateLibraryComplexity {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        input_file => {
            is  => 'String',
            doc => 'The SAM/BAM files to merge.  File type is determined by suffix.',
            is_many => 1,
        },
        output_file => {
            is  => 'String',
            doc => 'The output metrics file from picard..',
        },
        min_identical_bases => {
            doc => 'The minimum number of bases at the starts of reads that must be identical for reads to be grouped together for duplicate detection. In effect total_reads / 4^max_id_bases reads will be compared at a time, so lower numbers will produce more accurate results but consume exponentially more memory and CPU.',
            default_value => 5,
            is_optional => 1,
        },
        max_diff_rate => {
            doc => 'The maximum rate of differences between two reads to call them identical.',
            is_optional => 1,
            default_value => '0.03',
        },
        min_mean_quality => {
            doc => 'The minimum mean quality of the bases in a read pair for the read to be analyzed. Reads with lower average quality are filtered out and not considered in any calculations.',
            default_value => 20,
            is_optional => 1,
        },
        read_name_regex => {
            is => 'Text',
            doc => 'Regular expression that can be used to parse read names in the incoming SAM file. Read names are parsed to extract three variables: tile/region, x coordinate and y coordinate. These values are used to estimate the rate of optical duplication in order to give a more accurate estimated library size. The regular expression should contain three capture groups for the three variables, in order.',
            is_optional => 1,
            #default_value => '[a-zA-Z0-9\-\_]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*',
        },
        optical_duplicate_pixel_distance => {
            doc => 'The maximum offset between two duplicte clusters in order to consider them optical duplicates. This should usually be set to some fairly small number (e.g. 5-10 pixels) unless using later versions of the Illumina pipeline that multiply pixel values by 10, in which case 50-100 is more normal.',
            default_value => 100,
            is_optional => 1,
        },
    ],
};

sub help_brief {
    'Tool to estimate library complexity of raw reads(even unaligned) from a SAM/BAM file.';
}

sub help_detail {
    return <<EOS
    For Picard documentation of this command see:
    http://picard.sourceforge.net/command-line-overview.shtml#EstimateLibraryComplexity
EOS
}

sub execute {
    my $self = shift;

    my $complexity_cmd = $self->picard_path .'/EstimateLibraryComplexity.jar net.sf.picard.sam.EstimateLibraryComplexity';
    $complexity_cmd .= ' OUTPUT=' . $self->output_file;
    $complexity_cmd
        .= ' ' . join( ' ', map { "INPUT=$_" } $self->input_file );
    if ($self->min_identical_bases) {
        $complexity_cmd .= ' MIN_IDENTICAL_BASES='. $self->min_identical_bases;
    }
    if ($self->max_diff_rate) {
        $complexity_cmd .= ' MAX_DIFF_RATE='. $self->max_diff_rate;
    }
    if ($self->min_mean_quality) {
        $complexity_cmd .= ' MIN_MEAN_QUALITY='. $self->min_mean_quality;
    }
    if ($self->read_name_regex) {
        $complexity_cmd .= ' READ_NAME_REGEX='. $self->read_name_regex;
    }
    if ($self->optical_duplicate_pixel_distance) {
        $complexity_cmd .= ' OPTICAL_DUPLICATE_PIXEL_DISTANCE='. $self->optical_duplicate_pixel_distance;
    }
    $self->run_java_vm(
        cmd => $complexity_cmd,
        input_files => [$self->input_file],
        output_files => [$self->output_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
