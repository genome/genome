package Genome::Model::Tools::SvSim::FilterBreakdancer;

use Genome;
use strict;
use warnings;

my $BD_SCORE_FIELD = 8;

class Genome::Model::Tools::SvSim::FilterBreakdancer {
    is => "Command::V2",
    has_input => [
        input_file => {
            is => "Text",
            doc => "The input file containing breakdancer calls",
        },
        output_file => {
            is_output => 1,
            is => "Text",
            doc => "Path to write the filtered calls to",
        },
        min_score => {
            is => "Integer",
            doc => "Minimum bd score to pass filter",
        },
    ]
};

sub execute {
    my $self = shift;
    my $in = Genome::Sys->open_file_for_reading($self->input_file);
    my $out = Genome::Sys->open_file_for_writing($self->output_file);
    my $min_score = $self->min_score;

    while (my $line = $in->getline) {
        if ($line =~ /^#/) {
            $out->print($line);
            next;
        }

        chomp $line;
        my @fields = split("\t", $line);
        $out->print("$line\n") if $fields[$BD_SCORE_FIELD] >= $min_score;
    }
    return 1;
}

1;
