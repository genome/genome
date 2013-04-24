
package Genome::Model::Tools::Picard::ExtractIlluminaBarcodes;

# http://picard.sourceforge.net/command-line-overview.shtml#ExtractIlluminaBarcodes

use strict;
use warnings FATAL => 'all';

use Genome;

class Genome::Model::Tools::Picard::ExtractIlluminaBarcodes {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        basecalls_directory => {
            is  => 'String',
            doc => 'The basecalls output directory.',
        },
        output_directory => {
            is  => 'String',
            doc => 'Where to write _barcode.txt files.',
            is_optional => 1,
        },
        lane => {
            is  => 'String',
            doc => 'Lane number.',
        },
        read_structure => {
            is  => 'String',
            doc => 'A description of the logical structure of clusters in an Illumina Run',
        },
        barcode_file => {
            is  => 'String',
            doc => 'Tab-delimited file of barcode sequences, barcode name and and optionally library name.',
        },
        metrics_file => {
            is  => 'String',
            doc => 'Per-barcode and per-lane metrics written to this file.',
        },
        max_mismatches => {
            is  => 'String',
            doc => 'Maximum mismatches for a barcode to be considered a match.',
        },
        min_mismatch_delta => {
            is  => 'String',
            doc => 'Minimum difference between number of mismatches in the best and second best barcodes for a barcode to be considered a match.',
        },
        max_no_calls => {
            is  => 'String',
            doc => 'Maximum allowable number of no-calls in a barcode read before it is considered unmatchable.',
        },
        minimum_base_quality => {
            is  => 'String',
            doc => 'Minimum base quality. Any barcode bases falling below this quality will be considered a mismatch even in the bases match.',
            is_optional => 1,
        },
        compress_outputs => {
            is  => 'String',
            doc => 'Compress output s_l_t_barcode.txt files using gzip and append a .gz extension to the filenames.',
            is_optional => 1,
        },
        num_processors => {
            is  => 'String',
            doc => 'Run this many threads in parallel.',
            is_optional => 1,
        },
    ],
};

sub help_brief {
    'Determine the barcode for each read in an Illumina lane.',
}

sub help_detail {
    return <<EOS
    Determine the barcode for each read in an Illumina lane.  For Picard documentation of this command see:
    http://picard.sourceforge.net/command-line-overview.shtml#ExtractIlluminaBarcodes
EOS
}

sub execute {
    my $self = shift;

    my $jar_path = $self->picard_path . '/ExtractIlluminaBarcodes.jar';
    unless (-e $jar_path) {
        die('Failed to find '. $jar_path .'!  This command may not be available in version '. $self->use_version);
    }

    my %map_args = qw(
        basecalls_directory basecalls_dir
        output_directory    output_dir
    );

    my $args = join(' ',
        map {
            my $value = $self->$_;
            my $arg = $map_args{$_} || $_;
            defined($value) ? (uc($arg) . "='$value'") : ()
        } sort qw(
            basecalls_directory
            output_directory
            lane
            read_structure
            barcode_file
            metrics_file
            max_mismatches
            min_mismatch_delta
            max_no_calls
            minimum_base_quality
            compress_outputs
            num_processors
        )
    );

    my $cmd = $jar_path . " net.sf.picard.illumina.ExtractIlluminaBarcodes $args";
    $self->run_java_vm(
        cmd          => $cmd,
        input_files  => [ $self->barcode_file, ],
#        output_files => [ $self->metrics_file ],
#        skip_if_output_is_present => 0,
    );
    return 1;
}

1;
__END__

