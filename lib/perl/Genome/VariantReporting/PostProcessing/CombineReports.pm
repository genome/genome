package Genome::VariantReporting::PostProcessing::CombineReports;

use strict;
use warnings;
use Genome;
use List::Util qw(first);
use List::MoreUtils qw(firstidx);
use Set::Scalar;
use Memoize;

class Genome::VariantReporting::PostProcessing::CombineReports {
    is => 'Command',
    has => [
        reports => {
            is => 'Path',
            doc => 'The reports you wish to combine. They must all be the same type of report (same columns).',
            is_many => 1,
        },
        sort_columns => {
            is_optional => 1,
            is => 'Text',
            is_many => 1,
            doc => 'Column names to sort by. If the reports have no headers, provide column numbers (1-based). If not provided, the report will be unsorted.',
        },
        contains_header => {
            is => 'Boolean',
            doc => 'Set to true if the report contains headers'
        },
        output_file => {
            is => 'Path',
            doc => 'Where the combined report should go',
        },
        separator => {
            is => 'Text',
            default => "\t",
            doc => 'Field separator for the reports',
        },
        columns_to_split => {
            is => 'Text',
            is_optional => 1,
            is_many => 1,
            doc => 'Columns that should be split up.  These columns contain key:value pairs.  The key will be appended to the column header and the value will be put in the column.  Only valid if the file has a header.',
        },
    ],
};

sub execute {
    my $self = shift;

    $self->validate;

    my $combined_file = $self->combine_files;

    my $sorted_file = $self->sort_file($combined_file);

    if ($self->columns_to_split) {
        my $split_file = $self->split_file($sorted_file);
        $self->move_file_to_output($split_file);
    }
    else {
        $self->move_file_to_output($sorted_file);
    }
    return 1;
}

sub move_file_to_output {
    my $self = shift;
    my $file = shift;
    Genome::Sys->move_file(
        $file, $self->output_file
    );
}

sub combine_files {
    my $self = shift;

    my $combine_command;
    if ($self->contains_header) {
        $combine_command = 'tail -n +2 %s >> %s';
    } else {
        $combine_command = 'cat %s >> %s';
    }

    my $combined_file = Genome::Sys->create_temp_file_path;
    for my $report ($self->reports) {
        Genome::Sys->shellcmd(cmd => sprintf($combine_command, $report, $combined_file));
    }
    return $combined_file;
}

# Sort the file if required. Regardless, put the header in place.
sub sort_file {
    my ($self, $combined_file) = @_;

    my $sorted_file = Genome::Sys->create_temp_file_path;
    if (defined $self->sort_columns) {
        my $fh = Genome::Sys->open_file_for_writing($sorted_file);
        $self->print_header_to_fh($fh);
        Genome::Sys->shellcmd(cmd => sprintf('sort %s %s >> %s', $self->get_sort_params, $combined_file, $sorted_file));
    } else {
        my ($fh, $header_file) = Genome::Sys->create_temp_file;
        $self->print_header_to_fh($fh);
        Genome::Sys->concatenate_files( [$header_file, $combined_file], $sorted_file );
    }
    return $sorted_file;
}

sub split_file {
    my ($self, $file) = @_;
    my $split_file = Genome::Sys->create_temp_file_path;
    my $fh = Genome::Sys->open_file_for_writing($split_file);

    my %keys_to_append;
    my $in = Genome::Sys->open_file_for_reading($file);
    my $header = <$in>;
    chomp $header;
    my @header_fields = split($self->separator, $header);
    for my $header_field (@header_fields) {
        if (first {$_ eq $header_field} $self->columns_to_split) {
            $keys_to_append{$header_field} = {};
        }
    }
    while (my $line = <$in>) {
        chomp $line;
        my @fields = split ($self->separator, $line);
        for my $field_name (keys %keys_to_append) {
            my $field_index = firstidx {$_ eq $field_name} @header_fields; 
            my @values = split(",", $fields[$field_index]);
            for my $value (@values) {
                unless ($value eq ".") {
                    my ($sub_field, $sub_value) = split(":", $value);
                    $keys_to_append{$field_name}->{$sub_field} = $header_fields[$field_index]."-".$sub_field;
                }
            }
        }
    }
    $in->close;
    my @new_header;
    for my $header_field (@header_fields) {
        if (defined $keys_to_append{$header_field}) {
            for my $split_header (values %{$keys_to_append{$header_field}}) {
                push @new_header, $split_header;
            }
        }
        else {
            push @new_header, $header_field;
        }
    }
    $fh->print(join($self->separator, @new_header),"\n");
    $in = Genome::Sys->open_file_for_reading($file);
    $header = <$in>;
    while (my $line = <$in>) {
        chomp $line;
        my @fields = split ($self->separator, $line);
        my @new_fields;
        my $counter = 0;
        for my $header_field (@header_fields) {
            if (defined $keys_to_append{$header_field}) {
                my @split_field = split(",", $fields[$counter]);
                my %split_field_dict;
                for my $split_field_item (@split_field) {
                    my ($subfield, $subvalue) = split(":", $split_field_item);
                    $split_field_dict{"$header_field-$subfield"} = $subvalue;
                }
                for my $split_header (values %{$keys_to_append{$header_field}}) {
                    if (defined $split_field_dict{$split_header}) {
                        push @new_fields, $split_field_dict{$split_header};
                    }
                    else {
                        push @new_fields, ".";
                    }
                }
            }
            else {
                push @new_fields, $fields[$counter];
            }
            $counter++;
        }
        $fh->print(join($self->separator, @new_fields),"\n");
    }
    $fh->close;
    return $split_file;
}

sub print_header_to_fh {
    my ($self, $fh) = @_;
    if ($self->contains_header) {
        $fh->print( join($self->separator, $self->get_master_header) . "\n");
    }
    $fh->close;
}

# Make sure all inputs and outputs are readable. Make sure all headers are the same. Make sure sort_columns are contained in the header (this also ensures they are numeric if they must be).
sub validate {
    my $self = shift;

    if ($self->columns_to_split and !$self->contains_header) {
        die $self->error_message("If columns-to-split are specified, then a header must be present");
    }
    Genome::Sys->validate_file_for_writing($self->output_file);

    my $master_header = Set::Scalar->new($self->get_master_header);
    for my $report ($self->reports) {
        Genome::Sys->validate_file_for_reading($report);

        my $current_header = Set::Scalar->new($self->get_header($report));
        unless ($current_header->is_equal($master_header)) {
            die $self->error_message("Headers for the reports are not the same. First header:\n%s\nCurrent header:\n%s", $master_header, $current_header);
        }
    }

    my $sort_columns = Set::Scalar->new($self->sort_columns);
    unless($master_header->contains($sort_columns->members)) {
        die $self->error_message('The sort columns (%s) are not contained within the first header (%s)', $sort_columns, $master_header);
    }

    return 1;
}

sub get_sort_params {
    my $self = shift;
    return '-V ' . join " ", map { "-k$_" } $self->get_sort_column_numbers;
}

# Return the one-based indices of the columns by which we are sorting.
sub get_sort_column_numbers {
    my $self = shift;
    return unless (defined $self->sort_columns);

    # If the header is provided by name, we have to find the indices
    my @indices;
    if ($self->contains_header) {
        my @header = $self->get_master_header;
        for my $column ($self->sort_columns) {
            my $index = firstidx { $_ eq $column } @header;
            if ($index == -1) {
                die $self->error_message('Failed to find column (%s) in header (%s)', $column, join(",", @header));
            } else {
                push @indices, ($index+1); # We want 1-based numbers for sorting.
            }
        }
    } else {
        @indices = $self->sort_columns;
    }

    if (@indices) {
        return @indices;
    } else {
        die $self->error_message('Failed to get the indices for the sort columns (%s) in the master header (%s)', $self->sort_columns, join(",", $self->get_master_header) );
    }
}

sub get_master_header {
    my $self = shift;
    my @reports = $self->reports;
    return $self->get_header($reports[0]);
}
memoize('get_master_header');

# Given a file, return the header. If reports have no header, return an 'anonymous' one with just numbers.
sub get_header {
    my ($self, $file) = @_;

    my $fh = Genome::Sys->open_file_for_reading($file);
    my $line = $fh->getline;
    chomp $line;
    my @columns = split $self->separator, $line;

    if ($self->contains_header) {
        return @columns;
    } else {
        return ( 1..scalar(@columns) );
    }
}

1;
