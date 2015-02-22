package Genome::Model::Tools::SvSim::MergeBreakdancerResults;

use Genome;

use strict;
use warnings;

class Genome::Model::Tools::SvSim::MergeBreakdancerResults {
    is => "Command::V2",
    has_input => [
        input_files => {
            is => "Text",
            doc => "List of breakdancer input files to merge",
            is_many => 1,
        },

        output_file => {
            is => "Text",
            doc => "Merged output file",
            is_output => 1,
        },
    ]
};

sub execute {
    my $self = shift;

    my @inputs = map {Genome::Sys->open_file_for_reading($_)} $self->input_files;
    my $output = Genome::Sys->open_file_for_writing($self->output_file);

    for my $i (0..$#inputs) {
        while (my $line = $inputs[$i]->getline) {
            if ($line =~ /^#/) {
                $output->print($line) if $i == 0;
                next;
            }
            $output->print($line);
        }
        $inputs[$i]->close();
    }

    return 1;
}

1;
