package Genome::Model::SomaticValidation::Command::ValidateSvs;

use strict;
use warnings;

use Genome;

class Genome::Model::SomaticValidation::Command::ValidateSvs {
    is => 'Command::Tree',
    doc => 'validate SVs for somatic-validation models',
};

sub sub_command_category { 'pipeline steps' }

1;

