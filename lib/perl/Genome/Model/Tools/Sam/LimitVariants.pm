package Genome::Model::Tools::Sam::LimitVariants;

use strict;
use warnings;

use Genome;
use Genome::Model::Tools::RefCov::ROI::Bed;

class Genome::Model::Tools::Sam::LimitVariants {
    is => ['Command'],
    has => [
        variants_file => {
            is => 'Text',
            doc => 'A list of variants called by samtools',
        },
        bed_file => {
            is => 'Text',
            doc => 'A BED format file to limit variants within defined target regions',
        },
        output_file => {
            is => 'Text',
            doc => 'A file path to write the limited samtools variants to',
        }
    ],
    has_optional => [
        minimum_depth => {
            is => 'Integer',
            doc => 'The minimum depth to keep variant calls',
            default_value => 1,
        },
        wingspan => {
            is => 'Integer',
            doc => 'A wingspan to add to each target region to conserve variants within those expanded coordinates',
            default_value => 0,
        }
    ],
};

sub execute {
    my $self = shift;

    my $regions = Genome::Model::Tools::RefCov::ROI::Bed->create(
        file => $self->bed_file,
        region_index_substring => 5,
        wingspan => $self->wingspan,
        load_all => 1,
    );
    unless ($regions) {
        $self->error_message('Failed to create regions from BED file '. $self->bed_file);
        die($self->error_message);
    }

    my $output_fh = Genome::Sys->open_file_for_writing($self->output_file);
    unless ($output_fh) {
        $self->error_message('Failed to open output file for writing '. $self->output_file);
        die($self->error_message);
    }
    my $var_fh = Genome::Sys->open_file_for_reading($self->variants_file);
    unless ($var_fh) {
        $self->error_message('Failed to open variants file for reading '. $self->variants_file);
        die($self->error_message);
    }

    while (my $line = $var_fh->getline) {
        chomp($line);
        my @entry = split(/\s+/,$line);
        my $chrom = $entry[0];
        my $start = $entry[1];
        my $end = $entry[1];
        my $depth = $entry[7];
        if ($depth >= $self->minimum_depth && $regions->overlaps_regions($chrom,$start,$end)) {
            print $output_fh $line ."\n";
        }
    }
    $var_fh->close;
    $output_fh->close;
    return 1;
}

1;
