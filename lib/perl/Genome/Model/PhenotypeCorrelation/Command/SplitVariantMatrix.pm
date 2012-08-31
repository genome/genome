package Genome::Model::PhenotypeCorrelation::Command::SplitVariantMatrix;

use Digest::MD5;
use Genome;
use POSIX;
use Sort::Naturally qw/nsort/;

use strict;
use warnings;

class Genome::Model::PhenotypeCorrelation::Command::SplitVariantMatrix {
    is => "Genome::Command::Base",
    doc => "Split a variant matrix into multiple files by column",
    has => [
        input_file => {
            is => "File",
            doc => "The variant matrix to split",
        },
        output_prefix => {
            is => "File",
            doc => "Partial path for output files. The extension .#.txt will be appended to each output file where # is a number."
        },
        max_cols_per_file => {
            is => "Number",
            doc => "The maximum number of columns allowed in any given output file",
            default_value => 5000,
        },
    ],
};

sub execute {
    my ($self) = @_;

    my $delim = "\t";
    my $input_file = $self->input_file;
    my $prefix = $self->output_prefix;
    my $max_cols_per_file = $self->max_cols_per_file;

    die "Can't split an empty vcf matrix!" unless -s $input_file;
    my $ifh = Genome::Sys->open_file_for_reading($input_file);
    my $header_line = <$ifh>;
    chomp $header_line;
    my @header_fields = split($delim, $header_line);
    my $n_cols = scalar(@header_fields) - 1;
    # Maybe we only "split" into 1 file. Oh well
    my $n_files = int(ceil($n_cols / $max_cols_per_file));

    $self->status_message("Splitting $input_file: $n_cols columns into $n_files files...\n");

    # array of n_files arrayrefs containing [start, stop] column for each file
    # note that the column indices will start from 1 since we have row names
    my @file_ranges;
    my @ofhs; # output file handles
    my $out_paths = [];
    my $beg_col = 1;
    for (my $i = 0; $i < $n_files; ++$i) {
        my $out_path = "$prefix.$i.txt";
        push(@$out_paths, $out_path);
        my $ofh = Genome::Sys->open_file_for_writing("$out_path");
        push(@ofhs, $ofh);

        # note that end_col is the last col for file $i, not one past the last
        # to enable $beg_col..$end_col.
        my $end_col = $beg_col + $max_cols_per_file-1;
        $end_col = $n_cols if $end_col > $n_cols;
        push(@file_ranges, [$beg_col, $end_col]);

        $ofh->write(join("\t",
            $header_fields[0],
            @header_fields[$beg_col..$end_col]) . "\n");

        $beg_col = $end_col + 1;
    }

    while (<$ifh>) {
        chomp;
        my @fields = split($delim);
        # write row header to each file
        for (my $i = 0; $i < $n_files; ++$i) {
            # which columns does file $i want?
            my @columns = $file_ranges[$i]->[0]..$file_ranges[$i]->[1];
            $ofhs[$i]->write(
                join("\t",
                    $fields[0], # first write the row name
                    @fields[@columns] # then columns
                ) . "\n"
            );
        }
    }

    $ifh->close();

    return $out_paths;
}

1;
