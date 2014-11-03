package Genome::VariantReporting::Command::CombineReports;

use strict;
use warnings;
use Genome;
use List::Util qw(first);
use List::MoreUtils qw(firstidx);
use Set::Scalar;
use Memoize;

class Genome::VariantReporting::Command::CombineReports {
    is => 'Command',
    has_input => [
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
            is_output => 1,
        },
        separator => {
            is => 'Text',
            default => "\t",
            doc => 'Field separator for the reports',
        },
        split_indicators => {
            is => 'Text',
            is_optional => 1,
            is_many => 1,
            doc => 'A regular expression that indicates that columns whose headers match should be split up.  These columns contain key:value pairs.  The key will be appended to the column header and the value will be put in the column.  Only valid if the file has a header.',
        },
        entry_sources => {
            is => 'HASH',
            is_optional => 1,
            doc => 'Hash of report => TAG If entry_sources are specified, a column will be added to the combined report with a tag on each entry indicating which report it originally came from',
        },
    ],
};

sub execute {
    my $self = shift;

    my @reports_with_size = grep {-s $_} $self->reports;
    if (scalar(@reports_with_size) == 0) {
        #Create an empty output file
        Genome::Sys->touch($self->output_file);
        return 1;
    }

    $self->validate(@reports_with_size);

    my $combined_file = $self->combine_files(@reports_with_size);

    my $sorted_file = $self->sort_file($combined_file);

    if ($self->split_indicators) {
        my $split_file = $self->split_file($sorted_file);
        $self->move_file_to_output($split_file);
    }
    else {
        $self->move_file_to_output($sorted_file);
    }
    return 1;
}

sub columns_to_split {
    my $self = shift;
    return grep {
        my $field = $_;
        first {$field =~ $_} $self->split_indicators
    } $self->get_master_header;
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
    my @reports_with_size = @_;

    my $combined_file = Genome::Sys->create_temp_file_path;

    for my $report (@reports_with_size) {
        my $file_to_combine;
        if ($self->contains_header) {
            $file_to_combine = Genome::Sys->create_temp_file_path;
            my $reader = Genome::Utility::IO::SeparatedValueReader->create(
                input => $report,
                separator => $self->separator,
            );
            my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
                headers => [$self->get_master_header],
                print_headers => 0,
                separator => $self->separator,
                output => $file_to_combine,
            );
            while (my $entry = $reader->next) {
                $writer->write_one($entry);
            }
        } else {
            $file_to_combine = $report;
        }
        my $with_source = $self->add_source($report, $file_to_combine);
        my $combine_command = 'cat %s >> %s';
        Genome::Sys->shellcmd(cmd => sprintf($combine_command, $with_source, $combined_file));
    }
    return $combined_file;
}

sub add_source {
    my ($self, $report, $file) = @_;
    unless ($self->entry_sources) {
        return $file;
    }
    my $tag = $self->entry_sources->{$report};
    my $out_file = Genome::Sys->create_temp_file_path;
    my $out = Genome::Sys->open_file_for_writing($out_file);
    my $in = Genome::Sys->open_file_for_reading($file);

    while (my $line = <$in>) {
        chomp $line;
        print $out join($self->separator, $line, $tag)."\n";
    }
    $in->close;
    $out->close;
    return $out_file;
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

    my @header_fields = $self->get_header($file);
    my %keys_to_append = $self->get_keys_to_append($file);
    my @new_header = $self->calculate_new_header(\@header_fields, \%keys_to_append);

    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        headers => [@new_header],
        print_headers => 1,
        separator => $self->separator,
        output => $split_file,
    );

    my $reader = Genome::Utility::IO::SeparatedValueReader->create(
        separator => $self->separator,
        input => $file,
    );
    while (my $entry = $reader->next) {
        my $new_entry = $self->calculate_new_entry($entry, \@header_fields, \%keys_to_append);
        $writer->write_one($new_entry);
    }
    return $split_file;
}

sub get_keys_to_append {
    my ($self, $file) = @_;

    my %keys_to_append;
    for my $header_field ($self->columns_to_split) {
        $keys_to_append{$header_field} = {};
    }

    my $reader = Genome::Utility::IO::SeparatedValueReader->create(
        input => $file,
        separator => $self->separator,
    );

    while (my $entry = $reader->next) {
        for my $field_name ($self->columns_to_split) {
            my @values = split(",", $entry->{$field_name});
            for my $value (@values) {
                unless ($value eq ".") {
                    my ($sub_field, $sub_value) = split(":", $value);
                    $keys_to_append{$field_name}->{$sub_field} = $field_name."-".$sub_field;
                }
            }
        }
    }
    return %keys_to_append;
}

sub calculate_new_header {
    my ($self, $header_fields, $keys_to_append) = @_;
    my @new_header;
    for my $header_field (@$header_fields) {
        if (defined $keys_to_append->{$header_field}) {
            for my $split_header (values %{$keys_to_append->{$header_field}}) {
                push @new_header, $split_header;
            }
        }
        else {
            push @new_header, $header_field;
        }
    }
    return @new_header;
}

sub calculate_new_entry {
    my ($self, $entry, $header_fields, $keys_to_append) = @_;
    my $new_entry;
    for my $header_field (@$header_fields) {
        if (defined $keys_to_append->{$header_field}) {
            my @split_field = split(",", $entry->{$header_field});
            my %split_field_dict;
            for my $split_field_item (@split_field) {
                my ($subfield, $subvalue) = split(":", $split_field_item);
                $split_field_dict{"$header_field-$subfield"} = $subvalue;
            }
            for my $split_header (values %{$keys_to_append->{$header_field}}) {
                if (defined $split_field_dict{$split_header}) {
                    $new_entry->{$split_header} = $split_field_dict{$split_header};
                }
                else {
                    $new_entry->{$split_header} = ".";
                }
            }
        }
        else {
            $new_entry->{$header_field} = $entry->{$header_field};
        }
    }
    return $new_entry;
}

sub print_header_to_fh {
    my ($self, $fh) = @_;
    if ($self->contains_header) {
        $fh->print( join($self->separator, $self->get_master_header_with_source) . "\n");
    }
    $fh->close;
}

# Make sure all inputs and outputs are readable. Make sure all headers are the same. Make sure sort_columns are contained in the header (this also ensures they are numeric if they must be).
sub validate {
    my $self = shift;
    my @reports_with_size = @_;

    if ($self->split_indicators and !$self->contains_header) {
        die $self->error_message("If split_indicators are specified, then a header must be present");
    }

    if ($self->entry_sources) {
        for my $report ($self->reports) {
            my $tag = $self->entry_sources->{$report};
            unless (defined $tag) {
                die $self->error_message("No source tag defined for report $report");
            }
        }
    }

    Genome::Sys->validate_file_for_writing($self->output_file);

    my $master_header = Set::Scalar->new($self->get_master_header);
    for my $report (@reports_with_size) {
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

sub get_master_header_with_source {
    my $self = shift;
    if ($self->entry_sources) {
        return ($self->get_master_header, "Source");
    }
    else {
        return $self->get_master_header;
    }
}

# Given a file, return the header. If reports have no header, return an 'anonymous' one with just numbers.
sub get_header {
    my ($self, $file) = @_;

    my $reader = Genome::Utility::IO::SeparatedValueReader->create(
        input => $file,
        separator => $self->separator,
    );

    if ($self->contains_header) {
        return @{$reader->headers};
    } else {
        return ( 1..scalar(@{$reader->headers}) );
    }
}

1;
