package Genome::Model::Tools::BioSamtools::CompareAlignmentSummaries;

use strict;
use warnings;

use Genome;

my %sort_order = (
    label => 1,
    total_bp => 2,
    total_aligned_bp => 3,
    percent_aligned => 4,
    total_unaligned_bp => 5,
    total_duplicate_bp => 6,
    percent_duplicates => 7,
    paired_end_bp => 8,
    read_1_bp => 9,
    read_2_bp => 10,
    mapped_paired_end_bp => 11,
    proper_paired_end_bp => 12,
    singleton_bp => 13,
    total_target_aligned_bp => 14,
    percent_target_aligned => 15,
    unique_target_aligned_bp => 16,
    duplicate_target_aligned_bp => 17,
    percent_target_duplicates => 18,
    total_off_target_aligned_bp => 19,
    percent_off_target_aligned => 20,
    unique_off_target_aligned_bp => 21,
    duplicate_off_target_aligned_bp => 22,
    percent_off_target_duplicates => 23,
);

class Genome::Model::Tools::BioSamtools::CompareAlignmentSummaries {
    is => ['Genome::Model::Tools::BioSamtools'],
    has_input => [
        input_files => {
            is => 'Text',
            doc => 'A list of input alignment summaries'
        },
        output_file => {
            is => 'Text',
            doc => 'A file path to store tab delimited output',
        },
        labels => {
            doc => 'A list of text strings to be used as labels for each input file.  Must retain same order as input files list.',
            is_optional => 1,
        }
    ],
};

sub create {
    my $class = shift;
    my %params = @_;
    # retain sort order of input files and labels if provided
    if ($params{labels}) {
        my $labels = delete($params{labels});
        my $input_files = delete($params{input_files});
        my $self = $class->SUPER::create(%params);
        $self->input_files($input_files);
        $self->labels($labels);
        return $self;
    } else {
        return $class->SUPER::create(%params);
    }
}

sub execute {
    my $self = shift;

    my @data;
    my $i = 0;
    my @labels;
    if ($self->labels) {
        @labels = @{$self->labels};
    }
    my @headers;
    for my $input_file (@{$self->input_files}) {
        my $label;
        if ($self->labels) {
            $label = $labels[$i++];
        }
        my $reader = Genome::Utility::IO::SeparatedValueReader->create(
            separator => "\t",
            input => $input_file,
        );

        unless ($reader) {
            $self->error_message("Can't create SeparatedValueReader for input file $input_file");
            return;
        }
        while (my $data = $reader->next) {
            if ($label) {
                $data->{label} = $label;
            }
            push @data, $data;
            unless (@headers) {
                @headers = sort hash_sort_order (keys %{$data});
            }
        }
        $reader->input->close;
    }
    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        separator => "\t",
        headers => \@headers,
        output => $self->output_file,
    );
    for my $data (@data) {
        $writer->write_one($data);
    }
    $writer->output->close;
    return 1;
}

sub hash_sort_order {
    if (!defined($sort_order{$a}) || !defined($sort_order{$b})) {
        die('Failed to find sort order for column A '. $a .' or column B '. $b);
    }
    $sort_order{$a} <=> $sort_order{$b};
}

1;
