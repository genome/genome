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
        subclass_name => {
            is => 'String',
            is_optional => 1,
            valid_values => ['CollectRnaSeqMetrics', 'GenotypeConcordance', 'MarkDuplicates'],
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

    my $parser_class_name = 'Genome::Model::Tools::Picard';
    if (defined($self->subclass_name)) {
        $parser_class_name .= '::'. $self->subclass_name;
    }

    my @metrics_files = sort $self->metrics_files;

    my %metrics;
    for my $metrics_file (@metrics_files) {
         my $metrics_hash_ref = $parser_class_name->parse_file_into_metrics_hashref($metrics_file);
         if (defined($metrics{$metrics_file})) {
             $self->fatal_error_message('Duplicate metrics file: '. $metrics_file);
         }
         $metrics_hash_ref->{$self->additional_column_name} = $metrics_file;
         $metrics{$metrics_file} = $metrics_hash_ref;
    }

    my $first_metrics_hash_ref = $metrics{$metrics_files[0]};

    my @headers = keys %{$first_metrics_hash_ref};

    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $self->output_file,
        separator => $self->separator,
        headers => \@headers,
        ignore_extra_columns => $self->ignore_extra_columns,
    );
    unless ($writer) {
        $self->fatal_error_message('Unable to open output writer for file: '. $self->output_file);
    }
    for my $metrics_file (@metrics_files) {
        $writer->write_one($metrics{$metrics_file});
    }
    $writer->output->close();

    return 1;
}

1;
