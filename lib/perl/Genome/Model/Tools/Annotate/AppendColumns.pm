package Genome::Model::Tools::Annotate::AppendColumns;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Annotate::AppendColumns {
    is => "Command::V2",
    has => [
        input_variants => {
            is => 'File',
            doc => "File that contains the variants to annotate",
        },
        output_file => {
            is => 'File',
        },
        additional_columns_file => {
            is => 'File',
            doc => "File with the columns to append",
        },
        column_to_append => {
            is => 'Number',
            doc => "1-based index of the column to append.",
        },
        header => {
            is => 'Text',
            doc => "String to use as header for new column",
        },
        chrom_column => {
            is => 'Number',
            doc => '1-based index of the chromosome name column in additional-columns-file',
            default => 1,
        },
        start_column => {
            is => 'Number',
            doc => '1-based index of the start coord column in additional-columns-file',
            default => 2,
        },
        stop_column => {
            is => 'Number',
            doc => '1-based index of the stop coord column in the additional-columns-file',
            default => 3,
        },
    ],
};

sub execute {
    my $self = shift;

    my %additional_info;
    my %additional_info2;

    my $variants_in = Genome::Sys->open_file_for_reading($self->input_variants);
    my $variants_header_line = <$variants_in>;
    
    while(my $line = <$variants_in>) {
        chomp $line;
        my @fields = split(/\t/, $line);
        $additional_info2{$fields[0]}{$fields[1]}{$fields[2]} = 1;
    }
    $variants_in->close;

    my $in = Genome::Sys->open_file_for_reading($self->additional_columns_file);

    $self->status_message("Reading in additional columns");
    my $chrom_column = $self->chrom_column - 1;
    my $start_column = $self->start_column - 1;
    my $stop_column = $self->stop_column - 1;

    while(my $line = <$in>) {
        chomp $line;
        my @fields = split (/\t/, $line);
        if (defined $additional_info2{$fields[$chrom_column]}{$fields[$start_column]}{$fields[$stop_column]} and $additional_info2{$fields[$chrom_column]}{$fields[$start_column]}{$fields[$stop_column]} == 1) {
            $additional_info{$fields[$chrom_column]}{$fields[$start_column]}{$fields[$stop_column]} = $fields[$self->column_to_append - 1];
        }
    }
    $in->close;

    $in = Genome::Sys->open_file_for_reading($self->input_variants);
    my $out = Genome::Sys->open_file_for_writing($self->output_file);
    my $header_line = <$in>;
    chomp $header_line;
    $out->print(join("\t", $header_line, $self->header)."\n");

    $self->status_message("Writing final output file");
    while(my $line = <$in>) {
        chomp $line;
        my @fields = split (/\t/, $line);
        if (defined $additional_info{$fields[0]}{$fields[1]}{$fields[2]}) {
            $line = join ("\t", $line, $additional_info{$fields[0]}{$fields[1]}{$fields[2]});
        }
        else {
            $line = join ("\t", $line, "-");
        }
        $out->print($line."\n");
    }
    $in->close;
    $out->close;
}

1;

