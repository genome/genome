package Genome::Model::RnaSeq::Command; 

use strict;
use warnings;

use Genome;

class Genome::Model::RnaSeq::Command {
    is => 'Command::Tree',
    doc => 'operate on rna-seq models',
};

sub sub_command_category { 'type specific' }

1;
