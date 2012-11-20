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
        columns_to_append => {
            is => 'String',
            doc => "Headers of the columns to append.  Separate by commas.",
        },
    ],
};

sub execute {
    my $self = shift;

    my %additional_info;
    my $in = Genome::Sys->open_file_for_reading($self->additional_columns_file);

    my $header_line = <$in>;
    chomp $header_line;
    my @header_fields = split(/\t/, $header_line);

    my @columns_list = split ",", $self->columns_to_append;

    my %header;
    my $counter = 0;
    foreach my $header_field (@header_fields) {
        $header{$header_field} = $counter;
        $counter++;
    }
    while(my $line = <$in>) {
        chomp $line;
        my @fields = split (/\t/, $line);
        foreach my $column (@columns_list) {
            $additional_info{$fields[0]}{$fields[1]}{$fields[2]}{$column} = $fields[$header{$column}];
        }
    }
    $in->close;

    $in = Genome::Sys->open_file_for_reading($self->input_variants);
    my $out = Genome::Sys->open_file_for_writing($self->output_file);
    my $header_line = <$in>;
    chomp $header_line;
    $out->print(join("\t", $header_line, @columns_list)."\n");

    while(my $line = <$in>) {
        chomp $line;
        my @fields = split (/\t/, $line);
        foreach my $column (@columns_list) {
            if (defined $additional_info{$fields[0]}{$fields[1]}{$fields[2]}{$column}) {
                $line = join ("\t", $line, $additional_info{$fields[0]}{$fields[1]}{$fields[2]}{$column});
            }
            else {
                $line = join ("\t", $line, "-");
            }
        }
        $out->print($line."\n");
    }
    $in->close;
    $out->close;
}

1;

