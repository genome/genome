package Genome::Model::Tools::DetectVariants2::Result;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::DetectVariants2::Result {
    is => ['Genome::Model::Tools::DetectVariants2::Result::DetectionBase'],
    doc => 'This class represents the result of a variant detector.',
};

#This is a concrete subclass for "DetectionBase" that adds no additional functionality.

1;
