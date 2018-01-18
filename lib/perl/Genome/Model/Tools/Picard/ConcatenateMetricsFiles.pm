package Genome::Model::Tools::Picard::ConcatenateMetricsFiles;

use strict;
use warnings;

use Genome;
use Data::Dumper;

class Genome::Model::Tools::Picard::ConcatenateMetricsFiles {
    is  => 'Command::V2',
    has => [
        metrics_files   => {
            is  => 'String',
            doc => 'An input Picard metrics file',
            is_many => 1,
        },
        output_file  => {
            is  => 'String',
            doc => 'The output concatenated metrics file',
        },
        separator => {
            is => 'String',
            default_value => "\t",
        },
        ignore_extra_columns => {
            is => 'Boolean',
            default_value => 0,
        },
        additional_column_name => {
            is => 'String',
            default_value => 'SOURCE_FILE',
        },
    ],
};

sub help_brief {
    'Concatenate individual Picard metric files into a single file.';
}

sub help_detail {
    return <<EOS
    For Picard documentation of output metric files:
    https://broadinstitute.github.io/picard/picard-metric-definitions.html
EOS
}

sub execute {
    my $self = shift;

    my $writer;
    my %metrics_files;
    for my $metrics_file (sort $self->metrics_files) {
        if (defined($metrics_files{$metrics_file})) {
            $self->fatal_message('Duplicate metrics file %s', $metrics_file);
        }
        $metrics_files{$metrics_file} = 1;
        my $reader = Genome::Utility::IO::SeparatedValueReader->create(
            separator => "\t",
            ignore_lines_starting_with => '#|(?:^$)',
            input => $metrics_file,
        );
        $self->fatal_message('Failed to generate parser for %s', $metrics_file) unless $reader;
        unless ($writer) {
            my $headers = $reader->headers;
            my @writer_headers = @$headers;
            push @writer_headers, $self->additional_column_name;
            $writer = Genome::Utility::IO::SeparatedValueWriter->create(
                output => $self->output_file,
                separator => $self->separator,
                headers => \@writer_headers,
                ignore_extra_columns => $self->ignore_extra_columns,
            );
            unless ($writer) {
                $self->fatal_message('Unable to open output writer for %s', $self->output_file);
            }
        }
        while (my $metrics = $reader->next) {
            $metrics->{$self->additional_column_name} = $metrics_file;
            $writer->write_one($metrics);
        }
    }

    $writer->output->close();

    return 1;
}

1;
