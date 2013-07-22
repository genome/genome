package Genome::Model::Tools::Sword::Base;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sword::Base {
	is => 'Command::V2',
	is_abstract => 1,
};

sub help_detail {
    "This tree of tools will run Sword on RNA-seq BAM files.";
}

1;
