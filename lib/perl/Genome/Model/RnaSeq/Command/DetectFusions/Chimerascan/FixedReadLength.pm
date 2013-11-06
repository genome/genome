package Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::FixedReadLength;

use strict;
use warnings;

use Genome;

class Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::FixedReadLength {
    is => 'Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::Base',
};

sub command_class_prefix {
    return "Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::FixedReadLength";
}

