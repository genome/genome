package Genome::Model::Tools::BedTools::Intersect;

use strict;
use warnings;

use Genome;

my $DEFAULT_FILE_A_FORMAT = 'bam';
my $DEFAULT_FORCE_STRANDEDNESS = 0;
my $DEFAULT_INTERSECTION_TYPE = 'unique';

class Genome::Model::Tools::BedTools::Intersect {
    is => 'Genome::Model::Tools::BedTools',
    has_input => [
        input_file_a => {
            is => 'Text',
            doc => 'The file of lines with which to find overlaps.',
            shell_args_position => 1,
        },
        input_file_a_format => {
            is => 'Text',
            is_optional => 1,
            doc => 'The format of input file A',
            default_value => $DEFAULT_FILE_A_FORMAT,
            valid_values => ['bed','bam'],
        },
        input_file_b => {
            is => 'Text',
            doc => 'The file in which to find overlaps for each line of input-file-a.',
            shell_args_position => 2,
        },
        intersection_type => {
            is => 'Text',
            doc => 'The results to output: "a-only" returns those regions in file A but not overlapped in B; "unique" returns one line for each region in file A that is matched;"overlap_both" returns original A and B entries plus the number of base pairs of overlap between the two features. "v" means Only report those entries in A that have _no overlaps_ with B.',
            valid_values => ['a-only', 'unique', 'overlaps','overlap_both','write_both','v'],
            default_value => $DEFAULT_INTERSECTION_TYPE,
            is_optional => 1,
        },
        force_strandedness => {
            is => 'Boolean',
            is_optional => 1,
            default_value => $DEFAULT_FORCE_STRANDEDNESS,
            doc => 'Force strandedness.  That is, only include hits in A that overlap B on the same strand.'
        },
        minimum_overlap => {
            is => 'Number',
            is_optional => 1,
            doc => 'Minimum overlap required as a fraction of A.'
        },
        output_file => {
            is => 'Text',
            doc => 'The output file to write intersection output',
        },
        output_file_format => {
            is => 'Text',
            doc => 'The format to write the output file. "bam" only works if input_file_a is a BAM. Same as "input_file_a_format" by default.',
            valid_values => ['bed', 'bam'],
            is_optional => 1,
        },
        header => {
            is => 'Boolean',
            doc => 'Whether to reprint the header from file a',
            default_value => 0,
        },
    ],
};

sub help_brief {
    "Returns the depth and breadth of coverage of features from A on the intervals in B.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed-tools intersect --input-file-a example.bam --input-file-b other-list.bed --output-file intersection.bam --intersection-type all
EOS
}

sub help_detail {                           
    return <<EOS
More information about the BedTools suite of tools can be found at http://code.google.com/p/bedtools/. 
EOS
}

sub execute {
    my $self = shift;
    if ($self->header) {
        my $grep_cmd = 'grep "^#" '.$self->input_file_a." > ".$self->output_file;
        system($grep_cmd);
    }
    my $a_flag = '-a';
    if ($self->input_file_a_format eq 'bam') {
        $a_flag .= 'bam';
    }
    my $options = '';
    if ($self->force_strandedness) {
        $options .= ' -s';
    }
    if ($self->output_file_format) {
        if($self->input_file_a_format eq 'bed' and $self->output_file_format ne 'bed') {
            $self->error_message('Cannot output BAM from intersection of two BED files.');
            return;
        } elsif($self->input_file_a_format eq 'bam' and $self->output_file_format eq 'bed') {
            $options .= ' -bed';
        }
    }
    if (defined($self->minimum_overlap)) {
        $options .= ' -f '. $self->minimum_overlap;
    }
    if ($self->intersection_type) {
        if($self->intersection_type eq 'a-only') {
            $options .= ' -v';
        } elsif ($self->intersection_type eq 'unique') {
            $options .= ' -u';
        } elsif ($self->intersection_type eq 'overlap_both') {
            $options .= ' -wao';
        } elsif ($self->intersection_type eq 'write_both') {
            $options .= ' -wa -wb';
        } elsif ($self->intersection_type eq 'v') {
            $options .= ' -v';
        }
    }
    my $temp_file = Genome::Sys->create_temp_file_path;
    my $cmd = $self->bedtools_path .'/bin/intersectBed '. $options .' '. $a_flag .' '. $self->input_file_a .' -b '. $self->input_file_b .' > '. $temp_file;
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$self->input_file_a,$self->input_file_b],
        output_files => [$temp_file],
        skip_if_output_is_present => 0,
        allow_zero_size_output_files => 1,
    );
    my $cat_cmd = "cat $temp_file >> ".$self->output_file;
    system($cat_cmd);
    return 1;
}

1;
