package Genome::Model::Tools::Picard::SortSam;

use strict;
use warnings;

use Genome;
use IO::File;
use File::Basename;

my $DEFAULT_SORT_ORDER = 'coordinate';
my $DEFAULT_MAX_RECORDS_IN_RAM = 500000;


class Genome::Model::Tools::Picard::SortSam {
    is  => 'Genome::Model::Tools::Picard::Base',
    has_input => [
        input_file => {
            is  => 'String',
            doc => 'The SAM/BAM file to sort.  File type is determined by suffix.',
            picard_param_name => 'INPUT',
        },
        output_file => {
            is  => 'String',
            doc => 'The resulting sorted SAM/BAM file.  File type is determined by suffix.',
            picard_param_name => 'OUTPUT',
        },
        sort_order => {
            is => 'String',
            doc => 'The sort order of the merged output file.  default_value='. $DEFAULT_SORT_ORDER,
            valid_values => ['coordinate', 'unsorted', 'queryname'],
            default_value => $DEFAULT_SORT_ORDER,
            is_optional => 1,
            picard_param_name => 'SORT_ORDER',
        },
        max_records_in_ram => {
            doc => 'When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed.',
            is_optional => 1,
            default_value => $DEFAULT_MAX_RECORDS_IN_RAM,
            picard_param_name => 'MAX_RECORDS_IN_RAM',
        },
    ],
};

sub help_brief {
    'Tool to sort a BAM or SAM file using Picard';
}

sub help_detail {
    return <<EOS
    Tool to sort a BAM or SAM file using Picard.  For Picard documentation of this command see:
    http://broadinstitute.github.io/picard/command-line-overview.html#SortSam
EOS
}

sub _jar_name {
    return 'SortSam.jar';
}

sub _java_class_name {
    return 'net.sf.picard.sam.SortSam';
}

sub _shellcmd_extra_params {
    my $self = shift;
    return (
        input_files => [$self->input_file],
        output_files => [$self->output_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}

1;
