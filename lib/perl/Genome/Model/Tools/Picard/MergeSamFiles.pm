package Genome::Model::Tools::Picard::MergeSamFiles;

use strict;
use warnings;

use Genome;
use IO::File;
use File::Basename;

my $DEFAULT_SORT_ORDER = 'coordinate';
my $DEFAULT_ASSUME_SORTED = 0;
my $DEFAULT_MERGE_SEQUENCE_DICTIONARY = 1;

class Genome::Model::Tools::Picard::MergeSamFiles {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        input_files => {
            is  => 'List',
            doc => 'The SAM/BAM files to merge.  File type is determined by suffix.',
        },
        output_file => {
            is  => 'String',
            doc => 'The resulting merged SAM/BAM file.  File type is determined by suffix.',
        },
        assume_sorted => {
            is  => 'Integer',
            valid_values => [1, 0],
            doc => 'Assume the input data is sorted.  default_value='. $DEFAULT_ASSUME_SORTED,
            default_value => $DEFAULT_ASSUME_SORTED,
            is_optional => 1,
        },
        sort_order => {
            is => 'String',
            doc => 'The sort order of the merged output file.  default_value='. $DEFAULT_SORT_ORDER,
            valid_values => ['coordinate', 'unsorted', 'queryname'],
            default_value => $DEFAULT_SORT_ORDER,
            is_optional => 1,
        },
        merge_sequence_dictionary => {
            is => 'Integer',
            valid_values => [1, 0],
            doc => 'Merge the seqeunce dictionaries. Picard default is 0.  default_value='. $DEFAULT_MERGE_SEQUENCE_DICTIONARY,
            default_value => $DEFAULT_MERGE_SEQUENCE_DICTIONARY,
            is_optional => 1,
        },
    ],
};

sub help_brief {
    'Tool to merge BAM or SAM files using Picard';
}

sub help_detail {
    return <<EOS
    Tool to merge BAM or SAM files using Picard.  For Picard documentation of this command see:
    http://picard.sourceforge.net/command-line-overview.shtml#MergeSamFiles
EOS
}

sub execute {
    my $self = shift;
    my $input_files = $self->input_files;
    my @input_files;
    if (ref($input_files) eq 'ARRAY') {
        @input_files = @{$input_files};
    } else {
        @input_files = split(',',$input_files);
    }
    my $output_file = $self->output_file;

     my $merge_cmd = $self->picard_path .'/MergeSamFiles.jar net.sf.picard.sam.MergeSamFiles';
    if ($self->merge_sequence_dictionary) {
        $merge_cmd .= ' MSD=true';
    } else {
        $merge_cmd .= ' MSD=false';
    }
    $merge_cmd .= ' SO='. $self->sort_order;
    if ($self->assume_sorted) {
        $merge_cmd .= ' AS=true';
    } else {
        $merge_cmd .= ' AS=false';
    }
    $merge_cmd .= ' O='. $self->output_file;
    my $list_of_files = join(' I=',@input_files);
    $merge_cmd .= ' I='. $list_of_files;
    $self->run_java_vm(
        cmd => $merge_cmd,
        input_files => \@input_files,
        output_files => [$output_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
