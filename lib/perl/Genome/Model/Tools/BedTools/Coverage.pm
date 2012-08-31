package Genome::Model::Tools::BedTools::Coverage;

use strict;
use warnings;

use Genome;

my $DEFAULT_FILE_A_FORMAT = 'bam';
my $DEFAULT_HISTOGRAM = 0;
my $DEFAULT_FORCE_STRANDEDNESS = 0;

class Genome::Model::Tools::BedTools::Coverage {
    is => 'Genome::Model::Tools::BedTools',
    has_input => [
        input_file_a => {
            is => 'Text',
            doc => 'The A file when evaluating coverage, typically a BAM file.',
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
            doc => 'The B file when evaluating coverage.  Must be in BED format.',
        },
        histogram => {
            is => 'Boolean',
            is_optional => 1,
            default_value => $DEFAULT_HISTOGRAM,
            doc =>  'Report a histogram of coverage for each feature in B as well as a summary histogram for _all_ features in B',
        },
        force_strandedness => {
            is => 'Boolean',
            is_optional => 1,
            default_value => $DEFAULT_FORCE_STRANDEDNESS,
            doc => 'Force strandedness.  That is, only include hits in A that overlap B on the same strand.'

        },
        output_file => {
            is => 'Text',
            doc => 'The output file to write coverage output',
        },
    ],
};

sub help_brief {
    "Returns the depth and breadth of coverage of features from A on the intervals in B.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed-tools coverage ...
EOS
}

sub help_detail {                           
    return <<EOS
More information about the BedTools suite of tools can be found at http://code.google.com/p/bedtools/. 

Default Output:
         After each entry in B, reports:
           1) The number of features in A that overlapped the B interval.
           2) The number of bases in B that had non-zero coverage.
           3) The length of the entry in B.
           4) The fraction of bases in B that had non-zero coverage.

Histogram Output:
          Output (tab delimited) after each feature in B:
            1) depth
            2) # bases at depth
            3) size of B
            4) % of B at depth
    
EOS
}

sub execute {
    my $self = shift;

    print STDERR "disabled due to failing test, check and fix me...\n";
    return;

    my $a_flag = '-a';
    if ($self->input_file_a_format eq 'bam') {
        $a_flag .= 'bam';
    }
    my $options = '';
    if ($self->histogram) {
        $options .= ' -hist';
    }
    if ($self->force_strandedness) {
        $options .= ' -s';
    }
    my $cmd = $self->bedtools_path .'/bin/coverageBed '. $options .' '. $a_flag .' '. $self->input_file_a .' -b '. $self->input_file_b .' > '. $self->output_file;
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$self->input_file_a,$self->input_file_b],
        output_files => [$self->output_file],
    );
    return 1;
}

1;
