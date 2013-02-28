package Genome::Model::SomaticValidation::Command::ValidateLargeIndels;

use strict;
use warnings;

use Genome;

class Genome::Model::SomaticValidation::Command::ValidateLargeIndels {
    is => 'Command::Tree',
    doc => 'validate large indels for somatic-validation models',
};

sub sub_command_category { 'pipeline steps' }

1;

