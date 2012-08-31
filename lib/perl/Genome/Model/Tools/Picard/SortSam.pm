package Genome::Model::Tools::Picard::SortSam;

use strict;
use warnings;

use Genome;
use IO::File;
use File::Basename;

my $DEFAULT_SORT_ORDER = 'coordinate';
my $DEFAULT_MAX_RECORDS_IN_RAM = 500000;


class Genome::Model::Tools::Picard::SortSam {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        input_file => {
            is  => 'String',
            doc => 'The SAM/BAM file to sort.  File type is determined by suffix.',
        },
        output_file => {
            is  => 'String',
            doc => 'The resulting sorted SAM/BAM file.  File type is determined by suffix.',
        },
        sort_order => {
            is => 'String',
            doc => 'The sort order of the merged output file.  default_value='. $DEFAULT_SORT_ORDER,
            valid_values => ['coordinate', 'unsorted', 'queryname'],
            default_value => $DEFAULT_SORT_ORDER,
            is_optional => 1,
        },
        max_records_in_ram => {
            doc => 'When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed.',
            is_optional => 1,
            default_value => $DEFAULT_MAX_RECORDS_IN_RAM,
        },
    ],
};

sub help_brief {
    'Tool to sort a BAM or SAM file using Picard';
}

sub help_detail {
    return <<EOS
    Tool to sort a BAM or SAM file using Picard.  For Picard documentation of this command see:
    http://picard.sourceforge.net/command-line-overview.shtml#SortSam
EOS
}

sub execute {
    my $self = shift;

    my $jar_path = $self->picard_path .'/SortSam.jar';
    unless (-e $jar_path) {
        die('Failed to find '. $jar_path .'!  This command may not be available in version '. $self->use_version);
    }
    my $input_file = $self->input_file;
    my $output_file = $self->output_file;
    my $sort_cmd = $jar_path .' net.sf.picard.sam.SortSam O='. $output_file .' I='. $input_file .' SO='. $self->sort_order .' MAX_RECORDS_IN_RAM='. $self->max_records_in_ram;
    $self->run_java_vm(
        cmd => $sort_cmd,
        input_files => [$input_file],
        output_files => [$output_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
