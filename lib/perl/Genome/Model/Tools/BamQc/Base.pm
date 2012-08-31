package Genome::Model::Tools::BamQc::Base;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::BamQc::Base {
	is => 'Command::V2',
	is_abstract => 1,
};

sub help_detail {
    "These commands are intended to provide a robust analysis of BAM alignments combining many individual tools into a comprehensive evaluation.";
}

1;
