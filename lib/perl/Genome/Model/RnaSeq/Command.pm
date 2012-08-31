package Genome::Model::RnaSeq::Command; 

use strict;
use warnings;

use Genome;

class Genome::Model::RnaSeq::Command {
    is => 'Command::Tree',
};

sub sub_command_category { 'type specific' }

1;
