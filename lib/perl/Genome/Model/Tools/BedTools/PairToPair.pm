package Genome::Model::Tools::BedTools::PairToPair;

use strict;
use warnings;

use Genome;

my $DEFAULT_INTERSECTION_TYPE = 'both';

class Genome::Model::Tools::BedTools::PairToPair {
    is => 'Genome::Model::Tools::BedTools',
    has_input => [
        input_file_a => {
            is => 'Text',
            doc => 'The file of lines with which to find overlaps.',
            shell_args_position => 1,
        },
        input_file_b => {
            is => 'Text',
            doc => 'The file in which to find overlaps for each line of input-file-a.',
            shell_args_position => 2,
        },
        intersection_type => {
            is => 'Text',
            doc => "neither Report overlaps if neither end of A overlaps B. either  Report overlaps if either ends of A overlap B.  both    Report overlaps if both ends of A overlap B.  notboth Report overlaps if one or neither of A's overlap B",
            valid_values => ['neither', 'either', 'both', 'notboth'],
            default_value => $DEFAULT_INTERSECTION_TYPE,
            is_optional => 1,
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
        slop => {
            is => 'Integer',
            doc => 'The amount of slop (in b.p.). to be added to each footprint of A. *Note*: Slop is subtracted from start1 and start2 and added to end1 and end2.',
            default => 0,
        },
        slop_strand => {
            is => 'String',
            valid_values => ['+', '-', '+-'],
            doc => 'Add slop based to each BEDPE footprint based on strand.',
            default => '+-',
        },
        ignore_strand => {
            is => 'Boolean',
            doc => 'Ignore strands when searching for overlaps.',
            default => 0,
        },
        require_different_names => {
            is => 'Boolean',
            doc => 'Require the hits to have different names (i.e. avoid self-hits).',
            default => 0,
        },
    ],
};

sub help_brief {
    return "Report overlaps between two paired-end BED files (BEDPE).";
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed-tools pair-to-pair --input-file-a example.bedpe --input-file-b other-list.bedpe --output-file intersection.bedpe --intersection-type both
EOS
}

sub help_detail {                           
    return <<EOS
More information about the BedTools suite of tools can be found at http://code.google.com/p/bedtools/. 
EOS
}

sub execute {
    my $self = shift;
    my $options;

    if (defined $self->minimum_overlap) {
        $options .= " -f ".$self->minimum_overlap;
    }
    $options .= " -type ".$self->intersection_type;
    $options .= " -slop ".$self->slop;
    if ($self->slop_strand ne '+-') {
        $options .= " -ss ".$self->slop_strand;
    }
    if ($self->ignore_strand) {
        $options .= " -is";
    }
    if ($self->require_different_names) {
        $options .= " -rdn";
    }
    my $cmd = $self->bedtools_path .'/bin/pairToPair '. $options .' -a '.
        Genome::Sys->quote_for_shell($self->input_file_a) .
        ' -b '. Genome::Sys->quote_for_shell($self->input_file_b) .
        ' > '. Genome::Sys->quote_for_shell($self->output_file);
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$self->input_file_a,$self->input_file_b],
        output_files => [$self->output_file],
        skip_if_output_is_present => 0,
        allow_zero_size_output_files => 1,
    );

    return 1;
}

1;
