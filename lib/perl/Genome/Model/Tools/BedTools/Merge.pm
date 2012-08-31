package Genome::Model::Tools::BedTools::Merge;

use strict;
use warnings;

use Genome;

my $DEFAULT_MAX_DISTANCE = 1;
my $DEFAULT_REPORT_NUMBER = 0;
my $DEFAULT_REPORT_NAMES = 0;
my $DEFAULT_FORCE_STRANDEDNESS = 0;

class Genome::Model::Tools::BedTools::Merge {
    is => 'Genome::Model::Tools::BedTools',
    has_input => [
        input_file => {
            is => 'Text',
            doc => 'The input file to be merged',
        },
        output_file => {
            is => 'Text',
            doc => 'The output file to write coverage output',
        },
        maximum_distance => {
            is => 'Number',
            is_optional => 1,
            doc => 'Maximum distance between features allowed for features to be merged. For value 0(zero), overlapping & book-ended features are merged.
',
            default_value => $DEFAULT_MAX_DISTANCE,
        },
        report_number => {
            is => 'Boolean',
            is_optional => 1,
            default_value => $DEFAULT_REPORT_NUMBER,
            doc =>  'Report the number of BED entries that were merged. Note: "1" is reported if no merging occurred.',
        },
        report_names => {
            is => 'Boolean',
            is_optional => 1,
            default_value => $DEFAULT_REPORT_NAMES,
            doc => 'Report the names of the merged features separated by semicolons.',
        },
        force_strandedness => {
            is => 'Boolean',
            is_optional => 1,
            default_value => $DEFAULT_FORCE_STRANDEDNESS,
            doc => 'Force strandedness.  That is, only merge features that are the same strand.',
        },
    ],
};

sub help_brief {
    "Merges overlapping BED entries into a single interval.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed-tools merge ...
EOS
}

sub help_detail {
    return <<EOS
More information about the BedTools suite of tools can be found at http://code.google.com/p/bedtools/. 
EOS
}

sub execute {
    my $self = shift;

    my $options = ' -d '. $self->maximum_distance;
    if ($self->report_number) {
        $options .= ' -n';
    }
    if ($self->report_names) {
        $options .= ' -nms';
    }
    if ($self->force_strandedness) {
        $options .= ' -s';
    }
    my $cmd = $self->bedtools_path .'/bin/mergeBed '. $options .' -i '. $self->input_file .' > '. $self->output_file;
    # Let's turn off the status message printing to STDERR
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$self->input_file],
        output_files => [$self->output_file],
        skip_if_output_is_present => 0,
        print_status_to_stderr => 0,
    );
    return 1;
}

1;
