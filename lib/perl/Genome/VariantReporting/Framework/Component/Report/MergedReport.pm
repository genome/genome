package Genome::VariantReporting::Framework::Component::Report::MergedReport;

use strict;
use warnings;
use Genome;
use List::Util qw(first);
use List::MoreUtils qw(firstidx);
use Set::Scalar;
use Memoize;
use Params::Validate qw(validate_pos :types);
use JSON qw(from_json);

our $REPORT_PKG = 'Genome::VariantReporting::Framework::Component::Report::SingleFile';

class Genome::VariantReporting::Framework::Component::Report::MergedReport {
    is => [
        'Genome::VariantReporting::Framework::Component::Report::MergeCompatible',
    ],
    has_input => [
        report_results => {
            is => 'Genome::VariantReporting::Framework::Component::Report::MergeCompatible',
            is_many => 1,
        },
    ],
    has_param => [
        sort_columns => {
            is_optional => 1,
            is => 'Text',
            is_many => 1,
        },
        contains_header => {
            is => 'Boolean',
        },
        use_header_from => {
            is => 'Genome::VariantReporting::Framework::Component::Report::MergeCompatible',
            is_optional => 1,
        },
        separator => {
            is => 'Text',
        },
        split_indicators => {
            is => 'Text',
            is_optional => 1,
            is_many => 1,
        },
        entry_sources => {
            is_many => 'Text',
            is => 'Text',
            is_optional => 1,
        },
    ],
    has_transient_optional => [
        _master_header => {
            is => 'ARRAY',
        },
    ],
};

sub _run {
    my $self = shift;

    my @reports_with_size = $self->get_reports_with_size;
    if (scalar(@reports_with_size) == 0) {
        #Create an empty output file
        Genome::Sys->touch($self->_temp_output_file);
        return 1;
    }

    $self->validate(@reports_with_size);

    $self->merge_legend_files();

    my $merged_file = $self->merge_files(@reports_with_size);

    my $sorted_file = $self->sort_file($merged_file);

    if ($self->split_indicators) {
        my $split_file = $self->split_file($sorted_file);
        $self->move_file_to_output($split_file);
    }
    else {
        $self->move_file_to_output($sorted_file);
    }
    return 1;
}

sub get_reports_with_size {
    my $self = shift;

    my @reports_with_size;
    for my $result ($self->report_results) {
        if (-s $result->report_path) {
            push @reports_with_size, $result->report_path;
        }
    }
    return @reports_with_size;
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
        $file, $self->_temp_output_file
    );
}

sub merge_legend_files {
    my $self = shift;

    my @legend_paths = $self->_get_input_report_legend_paths();
    return unless(@legend_paths);

    my $headers;
    my $all_filters = Set::Scalar->new();
    for my $legend_path (@legend_paths) {
        ($headers, my $filters) = _get_headers_and_filters($legend_path);
        $all_filters = $all_filters + $filters;
    }

    my @lines = (
        'Headers',
        @$headers,
        'Filters',
        sort $all_filters->members(),
    );

    my $output_path = File::Spec->join($self->temp_staging_directory,
            $self->legend_file_name);
    my $fh = Genome::Sys->open_file_for_writing($output_path);
    $fh->write(join("\n", @lines));
    $fh->close();
    return $output_path;
}

sub _get_input_report_legend_paths {
    my $self = shift;

    my @legends;
    for my $result ($self->report_results) {
        if ($result->can('legend_path') &&
            defined($result->legend_path) &&
            -f $result->legend_path) {
            push @legends, $result->legend_path;
        }
    }
    return @legends;
}

sub _get_headers_and_filters {
    my $legend_path = shift;

    my $fh = Genome::Sys->open_file_for_reading($legend_path);
    my $first_line = $fh->getline();
    unless ($first_line =~ /^Headers/) {
        die sprintf("Legend file (%s) appears to be malformed", $legend_path);
    }

    my @headers;
    my $filters = Set::Scalar->new();
    my $mode = 'headers';
    while (my $line = $fh->getline()) {
        chomp($line);
        if ($line =~ /^Filters/) {
            $mode = 'filters';
            next;
        }

        if ($mode eq 'headers') {
            push @headers, $line;
        } elsif ($mode eq 'filters') {
            $filters->insert($line);
        }
    }
    return \@headers, $filters;
}

sub legend_file_name {
    my $self = shift;
    my @report_results = $self->report_results;
    my $first_result = $report_results[0];
    if ($first_result->can('legend_file_name')) {
        return $first_result->legend_file_name;
    } else {
        return;
    }
}

sub legend_path {
    my $self = shift;

    if (defined($self->legend_file_name)) {
        return File::Spec->join($self->output_dir, $self->legend_file_name);
    } else {
        return;
    }
}

sub merge_files {
    my $self = shift;
    my @reports_with_size = @_;

    my $merged_file = Genome::Sys->create_temp_file_path;

    for my $report (@reports_with_size) {
        my $file_to_merge;
        if ($self->contains_header) {
            $file_to_merge = Genome::Sys->create_temp_file_path;
            my $reader = Genome::Utility::IO::SeparatedValueReader->create(
                input => $report,
                separator => $self->separator,
            );
            my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
                headers => [$self->get_master_header],
                print_headers => 0,
                separator => $self->separator,
                output => $file_to_merge,
            );
            while (my $entry = $reader->next) {
                $writer->write_one($entry);
            }
        } else {
            $file_to_merge = $report;
        }
        my $with_source = $self->add_source($report, $file_to_merge);
        my $merge_command = 'cat %s >> %s';
        Genome::Sys->shellcmd(cmd => sprintf($merge_command, $with_source, $merged_file));
    }
    return $merged_file;
}

sub add_source {
    my ($self, $report, $file) = @_;
    unless ($self->has_entry_sources) {
        return $file;
    }
    my $tag = $self->get_entry_source($report);
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
    my ($self, $merged_file) = @_;

    my $sorted_file = Genome::Sys->create_temp_file_path;
    if ($self->has_sort_columns) {
        my $fh = Genome::Sys->open_file_for_writing($sorted_file);
        $self->print_header_to_fh($fh);
        Genome::Sys->shellcmd(cmd => sprintf('sort %s %s >> %s', $self->get_sort_params, $merged_file, $sorted_file));
    } else {
        my ($fh, $header_file) = Genome::Sys->create_temp_file;
        $self->print_header_to_fh($fh);
        Genome::Sys->concatenate_files( [$header_file, $merged_file], $sorted_file );
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

sub _temp_output_file {
    my $self = shift;
    return File::Spec->join($self->temp_staging_directory, $self->file_name);
}

sub report_path {
    my $self = shift;
    return File::Spec->join($self->output_dir, $self->file_name);
}

sub can_be_merged {
    return 1;
}

sub merge_parameters {
    my $self = shift;
    my @report_results = $self->report_results;
    my $result_class = $report_results[0]->class;
    return $result_class->merge_parameters;
}


sub file_name {
    my $self = shift;
    my @report_results = $self->report_results;
    return $report_results[0]->file_name;
}


# Make sure all inputs and outputs are readable. Make sure all headers are the same. Make sure sort_columns are contained in the header (this also ensures they are numeric if they must be).
sub validate {
    my $self = shift;
    my @reports_with_size = @_;

    if ($self->split_indicators and !$self->contains_header) {
        die $self->error_message("If split_indicators are specified, then a header must be present");
    }

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
    return unless ($self->has_sort_columns);

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

    unless (defined($self->_master_header)) {
        my $result;
        if (defined($self->use_header_from) &&
            -s $self->use_header_from->report_path) {
            $result = $self->use_header_from;
        } else {
            my @report_results = grep {-s $_->report_path}
                $self->report_results;
            if (@report_results) {
                $result = shift @report_results;
            } else {
                die "No files available to get_master_header";
            }
        }
        my @header = $self->get_header($result->report_path);
        $self->_master_header(\@header);
    }
    return @{$self->_master_header};
}

sub get_master_header_with_source {
    my $self = shift;
    if ($self->has_entry_sources) {
        return ($self->get_master_header, "Source");
    }
    else {
        return $self->get_master_header;
    }
}

# Given a file, return the header. If reports have no header, 
# return an 'anonymous' one with just numbers.
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

sub has_entry_sources {
    my $self = shift;
    my @entry_sources = $self->entry_sources;
    return scalar(@entry_sources);
}

sub has_sort_columns {
    my $self = shift;
    my @sort_columns = $self->sort_columns;
    return scalar(@sort_columns);
}

sub get_entry_source {
    my ($self, $report) = validate_pos(@_, 1, 1);

    for my $entry_source ($self->entry_sources) {
        my ($report_id, $tag) = split(/\|/, $entry_source);
        my $report_result = Genome::VariantReporting::Framework::Component::Report::MergeCompatible->get($report_id);

        if ($report_result->report_path eq $report) {
            return $tag;
        }
    }
    die sprintf("No entry source for report (%s)", $report);
}

sub category {
    my $self = shift;

    my @report_users = map { $_->users('label like' => 'report:%') } $self->report_results;
    my $category;
    for my $user (@report_users) {
        if ($user->label =~ /report:(.*)/) {
            my $metadata_json = $1;
            my $m = from_json($metadata_json);
            if (!defined($category)) {
                $category = $m->{category};
            }
            elsif ($category ne $m->{category}) {
                die $self->error_message("Categories of unmerged reports (%s) are not the same: (%s), (%s)", join(', ', map { $_->id } @report_users), $category, $m->{category});
            }
        }
    }

    return $category;
}

1;
