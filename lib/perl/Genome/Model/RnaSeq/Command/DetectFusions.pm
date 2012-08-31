package Genome::Model::RnaSeq::Command::DetectFusions;

use strict;
use warnings;

use Genome;

class Genome::Model::RnaSeq::Command::DetectFusions {
    is => 'Command::Tree',
    doc => 'fusion detectors'
};

sub sub_command_category { 'pipeline' }

sub sub_command_sort_position { 3 }

1;

