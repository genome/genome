package Genome::Model::SingleSampleGenotype::Command;

use strict;
use warnings;

use Genome;

class Genome::Model::SingleSampleGenotype::Command {
    is => 'Command::Tree',
    doc => 'operate on single-sample-genotype models and builds',
};

sub sub_command_category { 'type specific' }

1;

