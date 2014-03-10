package Genome::Model::Tools::EdgeR::FilterCountsFile;

use Genome;
use strict;
use warnings;

my $R_SCRIPT = __FILE__ . ".R";

class Genome::Model::Tools::EdgeR::FilterCountsFile {
    is => "Command::V2",
    has_input => [
        counts_file => {
            is => 'File',
        },
        counts_per_million => {
            is => 'Integer',
        },
        num_samples => {
            is => 'Integer',
        },
    ],
    has_output => [
        filtered_counts_file => {
            is => 'File',
        },
    ],
};

sub execute {
    my $self = shift;

    my $cmd = sprintf("Rscript %s --input-file %s --counts-per-million %d --num-samples %d --output-file %s",
        $R_SCRIPT,
        $self->counts_file,
        $self->counts_per_million,
        $self->num_samples,
        $self->filtered_counts_file,
    );
    
    return Genome::Sys->shellcmd(
        cmd => $cmd,
    );
};

1;
