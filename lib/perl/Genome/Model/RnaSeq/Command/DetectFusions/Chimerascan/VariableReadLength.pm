package Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::VariableReadLength;

use strict;
use warnings;

use Genome;

class Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::VariableReadLength {
    is => 'Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::Base',
};

sub command_class_prefix {
    return "Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::VariableReadLength";
}

