package Genome::Model::Tools::Sambamba::Sort;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sambamba::Sort {
    is => 'Genome::Model::Tools::Sambamba::Base',
    has_optional_param => [
        memory => {
            is => 'Number',
            doc => 'approximate total memory limit for all threads (by default 2GB)',
            tool_arg_name => 'memory-limit',
        },
        queryname_sort => {
            is => 'Boolean',
            doc => 'sort by read name instead of coordinate (lexicographical order)',
            tool_arg_name => 'sort-by-name',
        },
        natural_sort => {
            is => 'Boolean',
            doc => 'sort by read name instead of coordinate (so-called natural sort as in samtools)',
            tool_arg_name => 'natural-sort',
        },
    ],
    has_input => [
        input_file => {
            is => 'FilePath',
            doc => 'The input BAM file to sort.',
            tool_input_file => 1,
            tool_bare_arg_position => 1,
        },
        output_file => {
            is => 'FilePath',
            is_optional => 1,
            doc => 'output file name; if not provided, the result is written to a file with .sorted.bam extension',
            tool_arg_name => 'out',
        },
    ],
};

sub help_detail{
    ''
}

sub _tool_subcommand_name {
    return 'sort';
}

sub _resolve_output_files {
    my $self = shift;

    if ($self->output_file) {
        return $self->output_file;
    }
}

1;
