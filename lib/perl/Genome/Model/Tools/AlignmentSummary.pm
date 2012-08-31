package Genome::Model::Tools::AlignmentSummary;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::AlignmentSummary {
    is => ['Command::Tree'],
};

sub help_brief {
    "A wrapper to run a C++ tool called alignment_summary"
}

sub help_detail {
    "This tool generates a YAML text file of alignment summary metrics, including the amount of on/off target when providing a BED file with or without wingspan."
}

1;
