package Genome::VariantReporting::PostProcessing::CombineReports;

use strict;
use warnings;
use Genome;
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
    ],
};

sub execute {
    my $self = shift;

    $self->validate;

    my $combined_file = $self->combine_files;

    my $sorted_file = $self->sort_file($combined_file);

    return 1;
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
# TODO sort non-numeric chromosomes to the bottom. Could use sed to do this for just X and Y but other chromosomes require something else.
# TODO This is especially problematic since the sortable columns no not necessarily nclude chromosome.
# TODO This could be fixed by making sort_columns less flexible - could change to chrom_column and pos_column and be very specific and implement our own sort.
sub sort_file {
    my ($self, $combined_file) = @_;

    if (defined $self->sort_columns) {
        my $fh = Genome::Sys->open_file_for_writing($self->output_file);
        $self->print_header_to_fh($fh);
        Genome::Sys->shellcmd(cmd => sprintf('sort %s %s >> %s', $self->get_sort_params, $combined_file, $self->output_file));
    } else {
        my ($fh, $header_file) = Genome::Sys->create_temp_file;
        $self->print_header_to_fh($fh);
        Genome::Sys->concatenate_files( [$header_file, $combined_file], $self->output_file );
    }
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
    return '-n ' . join " ", map { "-k$_" } $self->get_sort_column_numbers;
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
